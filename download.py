import os
import urllib.request as ur
import util

from meta import Metadata


def download(meta: Metadata, retry: bool):
    """
    Download valid files.
    :param meta: contains the ENCODE metadata
    :param retry: retry downloading failed files as specified in the configuration file, if True
    """
    # output format for the terminal download status update
    d_format = '\t\tdownloaded: {0:>5,}/{1:,}, failed: {2:>5,}'                 # type: str
    # the number of valid downloads for the specified assay(s) and the number of completed downloads
    all_downloads = meta.get_download_file_number(meta.DOWNLOAD_VALID)          # type: int
    completed = meta.get_download_file_number(meta.DOWNLOAD_COMPLETE)           # type: int

    print('\talready downloaded: {0:>5,}/{1:,}'.format(completed, all_downloads))

    download_list = meta.get_download_list()                                    # type: list[tuple[str, str, str, str]]

    # if retry is not enabled, compute the previously failed file IDs
    if retry:
        failed_ids = set()                                                      # type: set[str]
    else:
        failed_ids = meta.get_failed_ids()                                      # type: set[str]

    for file_id, file_url, file_path, file_dir in download_list:
        # if retry is not enabled, skip the previously failed file IDs
        if file_id in failed_ids:
            continue

        # create the directory in case it does not exist yet
        os.makedirs(file_dir, exist_ok=True)

        # try to download the file and store the file accession if the download fails
        try:
            ur.urlretrieve(url=file_url, filename=file_path)
            completed += 1
        except Exception:
            failed_ids.add(file_id)
            continue

        print(d_format.format(completed, all_downloads, len(failed_ids)), end='\r')
    print(d_format.format(completed, all_downloads, len(failed_ids)))

    # increment the failed counter of the files that failed to download
    meta.add_failed_downloads(failed_ids)
    # update the download information in the metadata
    meta.update_download_status()
    # output the updated metadata
    meta.update_files()


def compute_invalid_downloads(cfg: dict, meta: Metadata, delete: bool):
    """
    Compute the files that are not valid downloads but were downloaded anyway.
    :param cfg: configuration dictionary
    :param meta: contains the ENCODE metadata
    :param delete: True if invalid files should be deleted, False otherwise
    """
    # files that should be/have been downloaded
    valid_downloads = meta.get_file_list(meta.DOWNLOAD_VALID)                       # type: set[str]
    # remove files that are valid but were not downloaded successfully
    valid_downloads.intersection(meta.get_file_list(meta.DOWNLOAD_COMPLETE))
    # all downloaded files
    downloaded_files = set()                                                        # type: set[str]

    # compute all downloaded files
    for root, _, files in os.walk(cfg['file paths']['data directory']):
        for file_name in files:
            downloaded_files.add(os.path.join(root, file_name))

    # compute the downloaded files that shouldn't have been downloaded
    invalid_downloads = downloaded_files.difference(valid_downloads)
    print('\tnumber of invalid downloads: {0:,}'.format(len(invalid_downloads)))

    if delete:
        # number of files to remove, number of successfully removed files
        n = len(invalid_downloads)                                                  # type: int
        i = 0                                                                       # type: int
        # file paths that could not be removed
        failed = set()                                                              # type: set[str]
        r_format = '\tremoved: {0:>5,}/{1:,}, failed: {2:>5,}'                      # type: str

        # try to remove the incorrectly downloaded files and update the terminal
        for file_path in invalid_downloads:
            if util.remove_file(file_path):
                i += 1
            else:
                failed.add(file_path)

            print(r_format.format(i, n, len(failed)), end='\r')
        print(r_format.format(i, n, len(failed)))

        # update the download status and the corresponding metadata files
        meta.update_download_status()
        meta.update_files()

    # only consider files that are valid and were successfully downloaded
    meta.filter_download()
