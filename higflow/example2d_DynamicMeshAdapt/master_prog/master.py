import time
import datetime
import shutil
import subprocess
from pathlib import Path
import psutil
import sys

# status enum
IDLE = 'idle'
RUNNING = 'running'
FINISHED = 'finished'
KILLED = 'killed'
WAITING_COOLDOWN = 'Cooldown in Effect. Waiting...'

# Configuration
MAX_INSTANCES = 10
SHARED_FOLDER = 'shared'
CHECK_INTERVAL = 1  # seconds
LAST_LOG_INTERVAL = 5  # seconds
FULL_LOG_INTERVAL = 30  # seconds
COOLDOWN_TIME = 60  # seconds
THRESHOLD_TIME_FOR_DEAD_INSTANCE = 10 * FULL_LOG_INTERVAL  # seconds


SUBMISSION_FOLDER = Path(SHARED_FOLDER) / 'submission'
RUNNING_FOLDER = Path(SHARED_FOLDER) / RUNNING
FINISHED_FOLDER = Path(SHARED_FOLDER) / FINISHED
# STOPPED_FOLDER = Path(SHARED_FOLDER) / 'stopped'
LOG_FOLDER = Path(SHARED_FOLDER) / 'logs'

# Ensure required directories exist
for folder in [SUBMISSION_FOLDER, RUNNING_FOLDER, FINISHED_FOLDER, LOG_FOLDER]:
    folder.mkdir(parents=True, exist_ok=True)

# gets instance id from the log files
instance_files = list(LOG_FOLDER.glob('instance_*.log'))
instance_ids = [int(file.stem.split('_')[1]) for file in instance_files] if len(instance_files)>0 else [0]
INSTANCE_ID = min(set(range(1, max(instance_ids) + 2)) - set(instance_ids))
print(f"Instance ID: {INSTANCE_ID}")

# Set up logging
full_log_file = LOG_FOLDER / f'instance_{INSTANCE_ID}.log'
last_log_file = LOG_FOLDER / f'instance_{INSTANCE_ID}_last.log'

def get_script_files():
    return list(SUBMISSION_FOLDER.glob('*.sh'))

def get_instance_files():
    return list(LOG_FOLDER.glob('instance_*_last.log'))

def log_and_write_to_last(message, forced=False):
    current_time_struct = time.localtime(time.time())
    current_time_asctime = time.asctime(current_time_struct)
    message_str = f'{current_time_asctime} - {message}\n'
    curr_time = time.time()

    if forced or (int(curr_time) % FULL_LOG_INTERVAL == 0):
        try:
            with open(full_log_file, 'a') as f:
                f.write(message_str)
        except IOError as e:
            print(f"IOError writing to full log file: {e}", file=sys.stderr)
    if forced or (int(curr_time) % LAST_LOG_INTERVAL == 0):
        try:
            with open(last_log_file, 'w') as f:
                f.write(message_str)
        except IOError as e:
            print(f"IOError writing to last log file: {e}", file=sys.stderr)

def log_status(status, script_name=None, start_time=None, forced=False):
    if script_name:
        current_time = time.time()
        duration = int(current_time - start_time) if start_time else 0
        hours, minutes, seconds = duration // 3600, (duration % 3600) // 60, duration % 60
        log_and_write_to_last(f"Status: {status}, Script: {script_name}, Duration: {hours}:{minutes}:{seconds}",
                              forced=forced)
    else:
        log_and_write_to_last(f"Status: {status}", 
                              forced=forced)

def get_alive_instance_filenames():
    instance_files = get_instance_files()
    alive_instances = []
    for file in instance_files:
        try:
            with open(file, 'r') as f:
                lines = f.readlines()
                if lines:
                    datetime_str = lines[-1].split(' - ')[0]
                    last_timestamp = datetime.datetime.strptime(datetime_str, '%a %b %d %H:%M:%S %Y')
                    time_difference = datetime.datetime.now() - last_timestamp
                    if time_difference.total_seconds() < THRESHOLD_TIME_FOR_DEAD_INSTANCE:
                        alive_instances.append(file)
        except IOError as e:
            print(f"IOError reading instance file {file}: {e}", file=sys.stderr)
        except (IndexError, ValueError) as e:
            print(f"Error parsing timestamp in file {file}: {e}", file=sys.stderr)
    return alive_instances


def is_this_checking_turn():
    instance_names = get_alive_instance_filenames()
    try:
        instance_ids = sorted([int(file.stem.split('_')[1]) for file in instance_names])
        position = instance_ids.index(INSTANCE_ID)
        curr_time = time.time()
        time_mod_N = int(curr_time) % (2*MAX_INSTANCES)
        return 2*position == time_mod_N
    except ValueError as e:
        print(f"ValueError in is_this_checking_turn: {e}", file=sys.stderr)
        return False


def get_free_instance_ids():
    instance_files = get_alive_instance_filenames()
    free_instances = []
    for file in instance_files:
        with open(file, 'r') as f:
            lines = f.readlines()
            if lines and IDLE in lines[-1]:
                instance_num = int(file.stem.split('_')[1])
                free_instances.append(instance_num)
    return free_instances

def execute_script(script_path):
    log_status(RUNNING, script_path.name, forced=True)
    print(f"{RUNNING} script: {script_path.name}")

    script_start_time = time.time()
    # check for name conflicts before mving
    scr_stem = script_path.stem
    script_name = script_path.name
    while (RUNNING_FOLDER / script_name).exists():
        scr_stem += '_n'
        script_name = f"{scr_stem}{script_path.suffix}"

    try:
        shutil.move(script_path, RUNNING_FOLDER / script_name)
    except (shutil.Error, OSError) as e:
        print(f"Error moving script to running folder: {e}", file=sys.stderr)
        log_status(IDLE, forced=True)
        return

    # Run the script
    try: 
        process = subprocess.Popen(['bash', RUNNING_FOLDER / script_name])
    except subprocess.SubprocessError as e:
        print(f"Error starting subprocess: {e}", file=sys.stderr)
        log_status(IDLE, forced=True)
        return
    
    while process.poll() is None:
        # check if script is still in the running folder - if not, kill it
        if not (RUNNING_FOLDER / script_name).exists():
            terminate_process_and_children(process)
            log_status(KILLED, script_name, script_start_time, forced=True)
            return
        time.sleep(CHECK_INTERVAL)
        log_status(RUNNING, script_name, script_start_time)
    process.wait()

    # Move script to finished folder
    try:
        shutil.move(RUNNING_FOLDER / script_name, FINISHED_FOLDER / script_name)
    except (shutil.Error, OSError) as e:
        print(f"Error moving script to finished folder: {e}", file=sys.stderr)
    log_status(FINISHED, script_name=script_name, start_time=script_start_time, forced=True)

# check for cooldown by checking for the duration of all running instances and if they are greater than the cooldown
def cooldown_ongoing():
    instance_files = get_alive_instance_filenames()
    for file in instance_files:
        with open(file, 'r') as f:
            lines = f.readlines()
            if lines and RUNNING in lines[-1]: # if file is running
                for info in lines[-1].split(','): # get duration according to logging format above
                    if 'Duration:' in info:
                        duration = info.split(': ')[1]
                        hours, minutes, seconds = map(int, duration.split(':'))
                        total_seconds = hours * 3600 + minutes * 60 + seconds
                        if total_seconds < COOLDOWN_TIME:
                            return True
    return False

def terminate_process_and_children(p):
    try:
        parent = psutil.Process(p.pid)
        for child in parent.children(recursive=True):
            child.kill()
        parent.kill()
    except psutil.NoSuchProcess as e:
        print(f"Process already terminated: {e}", file=sys.stderr)
    except psutil.AccessDenied as e:
        print(f"Access denied when terminating process: {e}", file=sys.stderr)

def main():
    log_status(IDLE, forced=True)
    while True:
        log_status(IDLE)
        scripts = get_script_files()
        if scripts:
            if cooldown_ongoing():
                log_status(WAITING_COOLDOWN)
                time.sleep(CHECK_INTERVAL)
                continue

            if is_this_checking_turn():
                if INSTANCE_ID in get_free_instance_ids(): # likely redundant
                    execute_script(scripts[0])
        time.sleep(CHECK_INTERVAL)

if __name__ == "__main__":
    main()
