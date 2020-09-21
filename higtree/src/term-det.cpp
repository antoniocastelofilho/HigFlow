#include <mpi.h>

extern "C" {
#include "term-det.h"
}

#include <cstring>
#include <cassert>
#include <climits>
#include <cstdint>
#include <vector>
#include <list>

using SendBuff = std::vector<uint8_t>;

struct PendingMsg
{
	SendBuff buff;
	int dest;
};

struct OutgoingMsg
{
	OutgoingMsg() = delete;
	OutgoingMsg(size_t size):
		buff(size)
	{}
	OutgoingMsg(SendBuff&& buff):
		buff(std::move(buff))
	{}

	SendBuff buff;
	MPI_Request req;
};

struct TermDetection
{
	MPI_Comm comm;
	unsigned weight;
	unsigned total_weight;
	std::vector<uint8_t> recv_buf;
	int rank;

	std::list<OutgoingMsg> out_msgs;
	std::list<PendingMsg> pending_msgs;
	bool weight_req_sent;

	~TermDetection() {
		assert(pending_msgs.empty());
		for(auto &sent: out_msgs) {
			MPI_Wait(&sent.req, MPI_STATUS_IGNORE);
		}
		MPI_Barrier(comm);
	}
};

static const int ctrl_tag = 2734234;
static const int user_tag = 6814683;

enum CtrlMsgType
{
	TERMINATION_ANNOUNCE,
	WEIGHT_RETURN,
	WEIGHT_REQUEST,
	WEIGHT_RESPONSE
};

TermDetection *
term_det_init(MPI_Comm comm)
{
	auto td = new TermDetection;
	td->comm = comm;

	int comm_size;
	MPI_Comm_size(comm, &comm_size);

	MPI_Comm_rank(comm, &td->rank);

	td->weight = UINT_MAX / comm_size;
	td->total_weight = comm_size * (UINT_MAX / comm_size);

	td->weight_req_sent = false;

	return td;
}

void
term_det_destroy(TermDetection *td)
{
	delete td;
}

static void
send_ctrl_msg(TermDetection *td, int dest, CtrlMsgType type, unsigned value)
{
	td->out_msgs.emplace_back(2 * sizeof(unsigned));
	auto &out = td->out_msgs.back();

	unsigned *out_buf = reinterpret_cast<unsigned*>(&out.buff[0]);
	out_buf[0] = type;
	out_buf[1] = value;

	MPI_Isend(out_buf, 2, MPI_UNSIGNED, dest, ctrl_tag, td->comm, &out.req);
}

static bool
try_termination(TermDetection *td)
{
	if(td->weight >= td->total_weight) {
		assert(td->weight == td->total_weight);

		int comm_size;
		MPI_Comm_size(td->comm, &comm_size);

		for(int i = 1; i < comm_size; ++i) {
			send_ctrl_msg(td, i, TERMINATION_ANNOUNCE, 0);
		}

		return true;
	}

	return false;
}

static void
clear_sent(TermDetection *td)
{
	auto i = td->out_msgs.begin();
	while (i != td->out_msgs.end())
	{
		int done;
		MPI_Test(&i->req, &done, MPI_STATUS_IGNORE);
		if (done) {
			i = td->out_msgs.erase(i);
		} else {
			++i;
		}
	}
}

static void
do_send_async(TermDetection *td, SendBuff&& msg, int dest, unsigned w)
{
	assert(td->weight >= w);
	td->weight -= w;

	// Assign weight to the initial bytes
	*reinterpret_cast<unsigned*>(&msg[0]) = w;

	// Add to sent list, so we can check later if the operation is done.
	td->out_msgs.push_back(std::move(msg));

	// The message has been moved to out_msgs, its final resting
	// place. Use the references from there to send.
	OutgoingMsg &out = td->out_msgs.back();

	// The actual send call
	MPI_Isend(&out.buff[0], out.buff.size(), MPI_BYTE, dest, user_tag, td->comm, &out.req);
}

static void
send_next_pending(TermDetection *td, unsigned w)
{
	assert(!td->pending_msgs.empty());

	auto &msg_data = td->pending_msgs.front();
	do_send_async(td, std::move(msg_data.buff), msg_data.dest, w);
	td->pending_msgs.pop_front();
}

static void
request_more_weight(TermDetection *td)
{
	assert(td->weight == 1);
	assert(!td->weight_req_sent);

	send_ctrl_msg(td, 0, WEIGHT_REQUEST, td->weight);
	td->weight = 0;
	td->weight_req_sent = true;
}

static void
try_send_all_pending(TermDetection *td)
{
	while(!td->pending_msgs.empty() && td->weight > 1) {
		unsigned w = td->weight / 2;

		send_next_pending(td, w);
	}

	if(td->rank != 0 && !td->pending_msgs.empty() && !td->weight_req_sent) {
		request_more_weight(td);
	}
}

class TerminateException
{};

static void *
single_recv_message(TermDetection *td, MPI_Status &s, int *from, int *msg_size)
{
	void *ret = nullptr;

	if(s.MPI_TAG == ctrl_tag) {
		CtrlMsgType type;
		unsigned weight;
		{
			unsigned buf[2];
			MPI_Recv(buf, 2, MPI_UNSIGNED, s.MPI_SOURCE,
				ctrl_tag, td->comm, MPI_STATUS_IGNORE);
			type = CtrlMsgType(buf[0]);
			weight = buf[1];
		}

		// Except for TERMINATION_ANNOUNCE, all other messages carries weight
		// as value, so we sum it now, for TERMINATION_ANNOUNCE value is zero.
		td->weight += weight;

		switch(type) {
		case WEIGHT_RETURN:
			// It is a weight return message.
			assert(td->rank == 0);
			assert(td->weight > weight); // No weight overflow
			assert(td->weight <= td->total_weight); // No weight created
			break;

		case TERMINATION_ANNOUNCE:
			// Terminated!
			assert(weight == 0);
			throw TerminateException();

		case WEIGHT_REQUEST:
			assert(td->rank == 0);
			assert(td->weight > 1);
			{
				unsigned w = td->weight / 2;
				send_ctrl_msg(td, s.MPI_SOURCE, WEIGHT_RESPONSE, w);
				td->weight -= w;
			}
			break;

		case WEIGHT_RESPONSE:
			td->weight_req_sent = false;
			break;

		default:
			assert(false);
		}

	} else {
		// Must be a user message...
		assert(s.MPI_TAG == user_tag);

		int buf_size;
		MPI_Get_count(&s, MPI_BYTE, &buf_size);

		// Find space to store the message
		if(unsigned(buf_size) > td->recv_buf.size()) {
			td->recv_buf.resize(buf_size);
		}

		*from = s.MPI_SOURCE;
		MPI_Recv(&td->recv_buf[0], buf_size, MPI_BYTE, *from, user_tag,
			td->comm, MPI_STATUS_IGNORE);
		unsigned *buf_as_unsigned = reinterpret_cast<unsigned*>(&td->recv_buf[0]);
		td->weight += *buf_as_unsigned;

		*msg_size = buf_size - sizeof(unsigned);
		ret = buf_as_unsigned + 1;
	}

	// The incoming message may have brought enough weight
	// to send the pending messages. Try to send them:
	try_send_all_pending(td);

	return ret;
}

void *
term_det_wait_for_msg(TermDetection *td, int *from, int *msg_size, bool may_terminate)
{
	// Test status of all sent messages.
	clear_sent(td);

	if(may_terminate) {
		// Before preparing for termination, what for processes with rank > 0
		// incurs in a new message, we non-blocking receive pending messages.
		for(;;) {
			int flag;
			MPI_Status s;
			MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, td->comm, &flag, &s);

			if(!flag) {
				// There are no pending message!
				break;
			}

			try {
				void *ret = single_recv_message(td, s, from, msg_size);

				if(ret) {
					return ret;
				}
			} catch(TerminateException) {
				// This must never happen here.
				assert(false);
				throw;
			}
		}
	} else if (td->rank > 0 && td->weight > 1) {
		// Send some weight back to ensure we will never deadlock.
		unsigned w = td->weight / 2;
		send_ctrl_msg(td, 0, WEIGHT_RETURN, w);
		td->weight -= w;
	}

	for(;;) {
		if(may_terminate && td->pending_msgs.size() <= td->weight) {
			// Send all remaining pending msgs, there is no
			// problem to be out of weight, because we may
			// terminate now.
			while(!td->pending_msgs.empty()) {
				send_next_pending(td, 1);
			}

			if(td->rank == 0) {
				if(try_termination(td))
					return nullptr;
			} else if(td->weight > 0 && !td->weight_req_sent) {
				// Build and send termination message,
				// with all local weight
				send_ctrl_msg(td, 0, WEIGHT_RETURN, td->weight);
				td->weight = 0;
			}
		}

		MPI_Status s;
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, td->comm, &s);

		try {
			void *ret = single_recv_message(td, s, from, msg_size);

			if(ret) {
				return ret;
			}
		} catch(TerminateException) {
			assert(may_terminate);

			// Asserts this process has the ammount of weight it
			// must have in order to terminate.
			assert(td->weight == 0);

			return nullptr;
		}
	}
}

void
term_det_send_asyncv(TermDetection *td, void const * const *bufs, const size_t *buf_sizes,
	int to)
{
	// Calculate the total size to send
	size_t size = 0;
	for(const size_t *s = buf_sizes; *s; ++s) {
		size += *s;
	}

	// Allocate space for all buffers + weight
	SendBuff out(size + sizeof(unsigned));

	// Concatenate the message in the rest of the bytes
	size = sizeof(unsigned);
	for(const size_t *s = buf_sizes; *s; ++s, ++bufs) {
		memcpy(&out[size], *bufs, *s);
		size += *s;
	}

	// Check if we have enough weight to send the message.
	// If not, add the message to the pending queue.
	if(td->weight > 1) {
		// Calculate weight to send:
		// 1/2 of what currently holding.
		unsigned w = td->weight / 2;

		do_send_async(td, std::move(out), to, w);
	} else {
		// We have only 1 weight left, can't send the message.
		// Add it to the list to be sent later.
		td->pending_msgs.push_back({std::move(out), to});

		if(td->rank != 0 && !td->weight_req_sent) {
			request_more_weight(td);
		}
	}
}

void
term_det_send_async(TermDetection *td, const void *buf, size_t size, int to)
{
	size_t buf_sizes[] = {size, 0};
	term_det_send_asyncv(td, &buf, buf_sizes, to);
}
