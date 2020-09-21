#ifndef HIGTREE_SERIALIZE
#define HIGTREE_SERIALIZE

#include "higtree.h"
#include "higtree-iterator.h"

//! Must have the same starting members in the same order as hig_cell
//! in order to memcpy to work.
typedef struct __attribute__((packed)) hig_serial_tree
{
	Point lowpoint;
	Point highpoint;
	int numcells[DIM];
} hig_serial_tree;

/*! Serializes a complete higtree in a array of hig_serial_tree.
 *
 * The @p buf variable must be an array of appropriate size,
 * i.e. its size must be at least the total number of nodes
 * in the tree.
 */
size_t hs_serialize(hig_serial_tree *buf, hig_cell *tree);

/*! Serializes a complete higtree in a array of hig_serial_tree, displacing its
 * postion.
 *
 * The @p buf variable must be an array of appropriate size,
 * i.e. its size must be at least the total number of nodes
 * in the tree.
 */
size_t hs_serialize_displaced(hig_serial_tree *buf, hig_cell *tree, CPPoint disp);

/*! Allocates and deserialize a complete higtree from
 * an array of hig_serial_tree.
 */
size_t hs_deserialize(hig_cell **tree, hig_serial_tree *buf);

/*! Stores the state of a multi-step serialization. */
typedef struct hig_serialization_ctx
{
	higcit_celliterator *it;
	Point displacement;
} hig_serialization_ctx;

/*! Starts a multi-step serialization of a tree.
 *
 * After calling this function, you may call hs_serialize_n()
 * successivelly to serialize the tree.
 *
 * Clean up with hs_serialization_finish() afterwards.
 */
void hs_serialization_start(hig_serialization_ctx *hss, hig_cell *tree);

/*! Starts a multi-step displaced serialization of a tree, at a displaced
 * position.
 *
 * After calling this function, you may call hs_serialize_n()
 * successivelly to serialize the tree. The serialized tree will be
 * displaced by @p disp.
 *
 * Clean up with hs_serialization_finish() afterwards.
 */
void hs_displaced_serialization_start(hig_serialization_ctx *hss,
	hig_cell *tree, CPPoint disp);

/*! Serializes at most @p count elements into hig_serial_tree.
 *
 * At every invocation of this function, at most @count of the remanining
 * nodes of the tree are serialized. Serialization context must have been
 * initialzied with hs_serialization_start().
 *
 * @returns the number of nodes serialized. This is @p count unless there
 * is no more nodes to serialize, in which case the value returned is shorter.
 * Zero is returned if the serialization is over.
 */
size_t hs_serialize_n(hig_serialization_ctx *hss,
	hig_serial_tree *buf, size_t count);

/*! Frees the resource allocated by hs_serialization_start(). */
void hs_serialization_finish(hig_serialization_ctx *hss);


/*! Stores the state of a multi-step deserialization. */
typedef struct hig_deserialization_ctx
{
	hig_cell *tree;
	higcit_celliterator *it;
} hig_deserialization_ctx;

/*! Starts a multi-step deserialization of a tree.
 *
 * After calling this function, you may call hs_deserialize_n()
 * successivelly to deserialize the tree.
 *
 * Clean up with hs_deserialization_finish() afterwards.
 */
void hs_deserialization_start(hig_deserialization_ctx *hsd);

/*! Deserializes at most @p count elements into hig_serial_tree.
 *
 * At every invocation of this function, at most @p count nodes are processed
 * from the array @p buf. Deserialization context must have been
 * initialzied with hs_deserialization_start().
 *
 * It is safe to pass a @p count greater than the number of nodes needed
 * to complete the tree, in this case, the extra elements will be ignored.
 *
 * @returns the number of nodes deserialized. This is @p count unless the
 * deserialization is complete, in which case the actualy used to deserialize
 * the tree.
 */
size_t hs_deserialize_n(hig_deserialization_ctx *hsd,
	hig_serial_tree *buf, size_t count);

/*! Returns the tree built from the deserialization. */
hig_cell* hs_deserialization_finish(hig_deserialization_ctx *hsd);

#endif
