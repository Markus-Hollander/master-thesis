from collections import defaultdict
from datetime import timedelta
import multiprocessing as mp
import numpy as np
import os
import shutil
import subprocess as sp
from time import time
from typing import Any, Callable
import yaml


def load_yaml(file_path: str) -> dict:
    """
    Loads a YAML file.
    """
    with open(file_path, 'r') as file:
        return yaml.load(stream=file, Loader=yaml.SafeLoader)


def write_yaml(file_path: str, yaml_dict: dict) -> None:
    """
    Writes the content of a dictionary into a YAML file.
    """
    with open(file_path, 'w') as file:
        yaml.dump(yaml_dict, stream=file)


def remove_empty_sub_dirs(directory: str) -> None:
    """
    Removes empty sub-directories of a given directory, as well as the directory itself if it becomes empty.
    :param directory: path to the root directory
    """
    # cannot remove sub-directories of a directory that does not exist
    if not os.path.isdir(directory):
        return

    # recursively remove empty sub-directories
    for content in sorted(os.listdir(directory)):
        remove_empty_sub_dirs(os.path.join(directory, content))

    # the directory was already deleted by removedirs in a recursive call
    if not os.path.isdir(directory):
        return

    # try to remove the directory only if it's empty
    if not os.listdir(directory):
        try:
            os.removedirs(directory)
        except IOError:
            pass


def remove_file(file_path: str) -> bool:
    """
    Remove the specified file.
    :param file_path: path to the to be removed file
    :return: True if the file was removed, False otherwise
    """
    # cannot remove a file that does not exist
    if not os.path.isfile(file_path):
        return False

    # try to remove the file and handle potential errors
    try:
        os.remove(file_path)
        return True
    except IOError:
        return False


def delete_create_dir(dir_path: str, delete: bool) -> None:
    """
    Makes sure the directory exists. Deletes the previous content if specified.
    :param dir_path: path of the directory to create/delete
    :param delete: True if the directory content should get deleted, False otherwise
    """
    if delete:
        shutil.rmtree(dir_path, ignore_errors=True)
    os.makedirs(dir_path, exist_ok=True)


def print_time(start, indent=0) -> None:
    """
    Prints the time passed.
    :param start: start time
    :param indent: number of tabs at the beginning of the line
    """
    print('{0}--> time taken: {1}'.format('\t' * indent, timedelta(seconds=round(time() - start))))


def ddl() -> dict[Any, list]:
    """
    Enables multiprocessing to pickle default dicts.
    """
    return defaultdict(list)


def dds() -> dict[Any, set]:
    """
    Enables multiprocessing to pickle default dicts.
    """
    return defaultdict(set)


def base_command(log_path: str, cmd_args: list[str], overwrite_log=False) -> bool:
    """
    Executes the given command.
    :param log_path: path to the log file
    :param cmd_args: list of command line arguments to intersect executed
    :param overwrite_log: start with an empty log if True
    :return: True if the command terminated without errors, False otherwise
    """
    # start with an empty log or add to an already existing log
    mode = 'w' if overwrite_log else 'a'                                                        # type: str

    with open(log_path, mode) as file:
        try:
            # execute the command, log the command line output
            sp.check_call(args=' '.join(cmd_args), stdout=file, stderr=file, shell=True)

            file.write('\nSUCCESS\n')
            return True
        # log the error if one occurs
        except Exception as e:
            file.write('{0}\n\n'.format(e))
            file.write('\nFAILURE\n')

            return False


def intersect(ref_path: str, result_bed: str, intersect_path: str, strand=False) -> None:
    """
    Runs bedtools intersect to compute the intersection of intervals in the results with intervals in the reference.
    :param ref_path: path to a file with reference intervals
    :param result_bed: path to a file with intervals in the results
    :param intersect_path: path to the output file with the intersecting intervals
    :param strand: True if the intervals are strand specific, False otherwise
    """
    with open(intersect_path, 'w') as out_file:
        if strand:
            args = ['bedtools', 'intersect', '-s', '-wo', '-a', ref_path, '-b', result_bed]
        else:
            args = ['bedtools', 'intersect', '-wo', '-a', ref_path, '-b', result_bed]
        sp.check_call(args=args, stdout=out_file, stderr=sp.PIPE)


def closest(ref_path: str, result_bed: str, closest_path: str) -> None:
    """
    Runs bedtools closest to compute the reference intervals that are the closest to the result intervals.
    :param ref_path: path to a file with reference intervals
    :param result_bed: path to a file with intervals in the results
    :param closest_path: path to the output file with the closest intervals
    """
    with open(closest_path, 'w') as out_file:
        args = ['bedtools', 'closest', '-io', '-t', 'all', '-d', '-a', result_bed, '-b', ref_path]
        sp.check_call(args=args, stdout=out_file, stderr=sp.PIPE)


def process_chromosomes(chrom: str) -> str:
    """
    Unifies the spelling of chromosome identifiers.
    """
    return chrom.upper().strip('CHR')


def process_strand(strand: str) -> str:
    """
    Unifies the spelling of strand indicators.
    """
    if strand in ['1', '+']:
        return '+'
    elif strand in ['-1', '-']:
        return '-'
    return strand


def get(dictionary: dict[str, Any], key: str, default: Any) -> Any:
    """
    Tries to retrieve the value corresponding to the given key. Casts the value to int if possible, or to a set. If the
    key does not exist in the dictionary or the value is empty, return the default vlaue.
    :param dictionary: dictionary from which to retrieve the specified value
    :param key: dictionary key for which to retrieve and process the value
    :param default: default value to return if the corresponding value could not be retrieved or processed
    :return: value for the given key
    """
    if key in dictionary:
        value = dictionary[key]

        if not value:
            return default

        try:
            return int(value)
        except ValueError:
            return set(value.split(';'))

    return default


def overlap(x_1: int, y_1: int, x_2: int, y_2: int) -> int:
    """
    Computes the length of the overlap between two intervals
    :param x_1: start of the 1st interval
    :param y_1: end of the 1st interval
    :param x_2: start of the 2nd interval
    :param y_2: end of the 2nd interval
    :return: length of the overlap
    """
    return max(0, min(y_1, y_2) - max(x_1, x_2) + 1)


def remove_zero_rows(table: list[list[float | int]]) -> list[list[float | int]]:
    """
    Removes all rows from a table that only contain zeros in all columns.
    :param table: list of table rows
    :return: the processed table
    """
    if not table:
        return table

    # iterate in reverse order so that del table[i] does not lead to incorrect deletions
    for i in reversed(range(len(table))):
        if all(x == 0 for x in table[i]):
            del table[i]

    return table


def remove_zero_columns(table: list[list[float | int]]) -> list[list[float | int]]:
    """
    Removes all columns from a table that only contains zeros in all rows.
    :param table: list of table rows, each row being a list of column entries
    :return: the processed table
    """
    if not table:
        return table

    # iterate in reverse order so that del table[i] does not lead to incorrect deletions
    for j in reversed(range(len(table[0]))):
        if all(table[i][j] == 0 for i in range(len(table))):
            for i in reversed(range(len(table))):
                del table[i][j]

    # iterate in reverse order so that del table[i] does not lead to incorrect deletions
    for i in reversed(range(len(table))):
        if not table[i]:
            del table[i]

    return table


def remove_zero(table: list[list[float | int]]) -> list[list[float | int]]:
    """
    Removes all rows and columns from a table that only contain zeros.
    :param table: list of table rows, each row being a list of column entries
    :return: the processed table
    """
    return remove_zero_columns(remove_zero_rows(table))


def has_zero(table: list[list[float | int]]) -> bool:
    """
    Checks if any row or column in the table only contains zeros.
    :param table: list of table rows, each row being a list of column entries
    :return: True if any row or column in the table only contains zeros, False otherwise
    """
    for row in table:
        if set(row) == {0}:
            return True

    for j in range(len(table[0])):
        if all(table[i][j] == 0 for i in range(len(table))):
            return True

    return False


def write_table(file_path: str, table: list[list[Any]], header_rows: int, number_mask=None) -> None:
    """
    Writes the table to a tab-separated file.
    :param file_path: file path to which to write the table
    :param table: list of table rows, each row being a list of column entries
    :param header_rows: number of header rows
    :param number_mask: indices of columns that contain larger numbers
    """
    number_mask = number_mask if number_mask else []                                                # type: list[int]

    fmt = '\t'.join(['{' + str(j) + ':,}' if j in number_mask else '{' + str(j) + '}' for j in range(len(table[0]))]) + '\n'

    with open(file_path, 'w') as file:
        for row_idx in range(len(table)):
            if row_idx < header_rows:
                file.write('\t'.join(table[row_idx]) + '\n')
            else:
                file.write(fmt.format(*table[row_idx]))


def p_value_summary(table: list[list[str | float]], content_start: int, p_val_columns=None):
    """
    Extracts the rows and columns from the given table that contain p-values and formats those p-values.
    :param table: list of table rows, each row being a list of column entries
    :param content_start: row index at which the p-value content starts in the table
    :param p_val_columns: indices of columns that contain p-values
    :return: table with formatted p-values
    """
    p_val_columns = p_val_columns if p_val_columns else []                                          # type: list[int]

    for i in range(content_start, len(table)):
        for j in p_val_columns:
            p = table[i][j]

            if np.isnan(p):
                p = 'n.a.'
            elif p < 0.001:
                p = '< 0.001'
            elif p < 0.01:
                p = '< 0.01'
            elif p < 0.05:
                p = '< 0.05'
            else:
                p = round(p, 2)

            table[i][j] = p

    return table


def multiprocess_tasks(func: Callable, tasks: list, message: str, indent: int):
    """
    Uses multiprocessing to apply the given function to the provided list of tasks, while timing the execution and
    printing the execution status.
    :param func: function to run with multiprocessing
    :param tasks: lists of tasks to pass to the function
    :param message: message to print when the function starts and ends
    :param indent: number of tabs in front of the message
    :return: list of function call results
    """
    if not tasks:
        return []

    message = message.format(len(tasks)) if '{' in message else message

    start = time()
    print('\t' * indent + message, end='\r')
    with mp.Pool(maxtasksperchild=1) as pool:
        if isinstance(tasks[0], tuple):
            results = pool.starmap(func, tasks)
        else:
            results = pool.map(func, tasks)

    print('\t' * indent + message + ': done')
    print_time(start, indent=indent)
    return results


def execute(func: Callable, args, message: str, indent: int, override=True, removed=-1):
    """
    Executes the function with the given arguments, times the execution time and prints the execution status.
    :param func: the function to execute
    :param args: the arguments to pass to the function
    :param message: message to print when the function starts and ends
    :param indent: number of tabs in front of the message
    :param override: True if the start message should be replaced by the end message once the function is finished,
                     False if the end message should be printed separately
    :param removed: index of removed elements in the results
    :return: results of the function
    """
    message = message.format(len(args)) if '{' in message and args and not isinstance(args, tuple) else message

    start = time()
    if override:
        print('\t' * indent + message, end='\r')
    else:
        print('\t' * indent + message)

    if not args:
        results = func()
    elif isinstance(args, tuple):
        results = func(*args)
    else:
        results = func(args)

    if removed > -1:
        print('\t' * indent + message + ': done, removed: {0}'.format(len(results[removed])))
    else:
        print('\t' * indent + message + ': done')

    print_time(start, indent=indent)

    return results
