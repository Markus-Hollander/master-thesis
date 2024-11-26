import os
import pandas as pd
import sqlite3


class Metadata:
    DOWNLOAD_INCORRECT = 'incorrectly downloaded'
    DOWNLOAD_LIST = 'to download'
    DOWNLOAD_COMPLETE = 'download completed'
    DOWNLOAD_VALID = 'valid download'
    DOWNLOAD_FAILED = 'download fails'
    DOWNLOAD_REASON = 'download invalid reason'

    def __init__(self, in_path: str, download_config: dict, name_components=None, download_dir=None, filter_dict=None,
                 adjustment_dict=None):
        """
        :param in_path: path to the metadata .tsv file
        :param download_config: dictionary with information on what kind of files are valid to download
        :param name_components: list of _build_data_frame to include in the file name
        :param download_dir: path to the download directory
        :param filter_dict: dictionary with information on what file entries to include or exclude
        :param adjustment_dict: dictionary with information on how to adjust certain column entries
        """
        # construct a data frame from the metadata input .tsv
        self.df = pd.read_csv(in_path, dtype=str, na_values=[], keep_default_na=False, sep='\t')    # type: pd.DataFrame
        # adjust the column names by casting them to all lower case letters
        self.df.rename(columns=lambda c: c.lower().replace('_', ' '), inplace=True)
        # rename 2 frequently used _build_data_frame
        self.df.rename(columns={'biosample term name': 'cell', 'biosample organism': 'organism'}, inplace=True)
        
        print(self.df)

        # adjust the file entries if specifications are given
        if adjustment_dict:
            self._adjust(adjustment_dict)

        # filter the metadata if specifications are given
        if filter_dict:
            self._filter(filter_dict)

        # add or update file paths if specifications are given
        #if download_dir and name_components:
        #    self._add_file_paths(download_dir, name_components)

        # adjust the data type of numeric columns
        self._initialise_int()
        # initialise or update columns for downloading files
        self._initialise_download(download_config)

        self.db_path = ''
        self.tsv_path = ''

    def _initialise_int(self):
        """
        Replace empty strings with zeros and cast the column type to integer for the specified columns.
        """
        for column in ['size', 'read length']:
            self.df[column] = self.df[column].replace({'': 0}).astype(int)

    def _filter(self, filter_dict: dict[str, dict[str, dict[str, list[str]]]]):
        """
        Filter the metadata.
        """
        for filter_type, filters in filter_dict.items():
            for column, values in filters.items():
                # remove entries that do not contain at least one of the valid contents for specific columns
                if filter_type == 'include':
                    self.df = self.df[self.df[column].isin(values)]
                # remove entries that contain unwanted column entries
                elif filter_type == 'exclude':
                    self.df = self.df[~self.df[column].isin(values)]

    def _adjust(self, adjustment_dict: dict[str, dict]):
        """
        Adjust file entries.
        :param adjustment_dict: dictionary with information on how to adjust certain column entries
        """
        # adjust cell/tissue names
        for old, new in [(',', ''), (' ', '.'), ('-', '.'), ('_', '.')]:
            self.df['cell'] = self.df['cell'].str.replace(old, new)

        # adjust the run type column
        for old, new in [('single-ended', 'single'), ('paired-ended', 'paired')]:
            self.df['run type'] = self.df['run type'].str.replace(old, new)

        for filter_type, filters in adjustment_dict.items():
            for column, adjustments in filters.items():
                # remove unwanted sub-strings from the specified column entries
                if filter_type == 'remove':
                    for old in adjustments:
                        self.df[column] = self.df[column].str.replace(old, '')
                # replace entries with alternatives in the specified column entries
                elif filter_type == 'mapping':
                    self.df[column] = self.df[column].replace(adjustments)

    def _add_file_paths(self, download_dir: str, name_components: list[str]):
        """
        Adds or updates download directory, file name and file path information for each file.
        :param download_dir: path to the download directory
        :param name_components: list of columns that should go into the file name
        """
        def compute_dir(row: pd.Series) -> str:
            """
            :param row: a row in the metadata data frame
            :return: directory path
            """
            print(row['assay'])
            if row['assay'] == 'ChIP.seq':
                return os.path.join(download_dir, row['experiment target'], row['cell'])
            return os.path.join(download_dir, row['assay'], row['cell'])

        def compute_file_name(row: pd.Series) -> str:
            """
            :param row: a row in the metadata data frame
            :return: file name
            """
            file_format = row['file format']
            if len(file_format.split()) > 1:
                one, two = file_format.split()
                two = two.replace('bed', '').capitalize()
                file_format = one + two

            components = [str(row[col]) for col in name_components]

            return '_'.join(components).strip('_') + '.' + file_format

        self.df['file directory'] = self.df.apply(lambda row: compute_dir(row), axis=1)
        self.df['file name'] = self.df.apply(lambda row: compute_file_name(row), axis=1)
        self.df['file path'] = self.df.apply(lambda row: os.path.join(row['file directory'], row['file name']), axis=1)

    def _download_complete(self):
        """
        Computes for each file if it was downloaded completely.
        """
        def check_complete(row: pd.Series) -> bool:
            """
            :param row: a row in the metadata data frame
            :return: True if the file was downloaded and has the same size as specified in the metadata file
            """
            try:
                size = os.path.getsize(row['file path'])
                return size == row['size']
            except OSError:
                return False

        self.df[self.DOWNLOAD_COMPLETE] = self.df.apply(lambda row: check_complete(row), axis=1)

    def _valid_download(self, config: dict):
        """
        Adds or updates a column that contains True if a file is valid to download.
        :param config: dictionary with specifications of what constitutes a valid file
        """
        # extract the different specifications from the configuration dictionary
        disallowed_warnings = config['disallowed audit warnings']
        disallowed_internal = config['disallowed audit internal action']
        allowed_status = config['allowed file status']
        allowed_errors = config['allowed audit errors']
        allowed_not_compliant = config['allowed audit not compliant']
        allowed_formats = config['allowed formats']

        def check_valid(row: pd.Series) -> str:
            """
            :param row: a row representing a file in the metadata data frame
            :return: empty string if the file passes all requirements; reason for failure otherwise
            """
            if ',' in row['paired end']:
                return ', in paired end'

            if row['file status'] not in allowed_status:
                return 'not released'

            if set(row['audit warning'].split(', ')).intersection(disallowed_warnings):
                return 'disallowed warning'

            if set(row['audit internal action'].split(', ')).intersection(disallowed_internal):
                return 'disallowed internal action'

            if set(row['audit error'].split(', ')).difference(allowed_errors):
                return 'disallowed audit error'

            if set(row['audit not compliant'].split(', ')).difference(allowed_not_compliant):
                return 'disallowed audit not compliant'

            if row['assay'] not in allowed_formats.keys():
                return 'disallowed assay'

            if row['output type'] not in allowed_formats[row['assay']].keys():
                return 'disallowed output type'

            if row['file format'] not in allowed_formats[row['assay']][row['output type']]:
                return 'disallowed file format'

            return ''

        self.df[self.DOWNLOAD_REASON] = self.df.apply(lambda row: check_valid(row), axis=1)
        self.df[self.DOWNLOAD_VALID] = self.df.apply(lambda row: not row[self.DOWNLOAD_REASON], axis=1)

    def _download_list(self):
        """
        Adds or updates a column that contains True if a file is valid for download but has not been downloaded yet.
        """
        d_list, complete, valid = (self.DOWNLOAD_LIST, self.DOWNLOAD_COMPLETE, self.DOWNLOAD_VALID)
        self.df[d_list] = self.df.apply(lambda row: row[valid] and not row[complete], axis=1)

    def _invalid_downloads(self):
        """
        Adds or updates a column that contains True if a file has been downloaded completely but is not a valid file
        to download.
        """
        incorrect, complete, valid = (self.DOWNLOAD_INCORRECT, self.DOWNLOAD_COMPLETE, self.DOWNLOAD_VALID)
        self.df[incorrect] = self.df.apply(lambda row: row[complete] and not row[valid], axis=1)

    def _initialise_download(self, config: dict):
        """
        sets up and initialises columns relevant for file downloads.
        :param config: specifies what kind of files are valid to download
        """
        # add a column containing the number of failed download attempts, if it does not exist already
        if self.DOWNLOAD_FAILED not in self.df.keys():
            self.df[self.DOWNLOAD_FAILED] = 0
        self.df[self.DOWNLOAD_FAILED] = self.df[self.DOWNLOAD_FAILED].astype(int)

        # compute which files are valid to download
        self._valid_download(config)
        # compute which files have already been downloaded
        #self._download_complete()
        # compute which files have already been downloaded even though they are not valid
        #self._invalid_downloads()
        # compute which files are valid and still have to be downloaded
        #self._download_list()

    def update_download_status(self):
        """
        Recomputes which files have already been downloaded and which still need to be.
        """
        self._download_complete()
        self._download_list()

    def _get_filtered_download_frame(self, key: str) -> pd.DataFrame:
        """
        Filters the metadata down to entries where the specified download column is True.
        :param key: name of one of the boolean download columns
        :return: filtered metadata data frame
        """
        if key not in [self.DOWNLOAD_VALID, self.DOWNLOAD_COMPLETE, self.DOWNLOAD_INCORRECT, self.DOWNLOAD_LIST,
                       self.DOWNLOAD_FAILED]:
            raise ValueError('The column must be one of the boolean download columns.')

        # return the filtered metadata data frame
        if key != self.DOWNLOAD_FAILED:
            return self.df[self.df[key]]

        # the failed download column contains integers, not booleans
        return self.df[self.df[key] > 0]

    def get_download_file_number(self, key: str) -> int:
        """
        Computes the number of valid files where the specified download column is True.
        :param key: name of one of the boolean download columns
        :return: number of files
        """
        filtered = self._get_filtered_download_frame(key)

        # compute the number of files that are valid but were not completely downloaded
        if key == self.DOWNLOAD_FAILED:
            return len(filtered[filtered[self.DOWNLOAD_VALID]
                                & ~filtered[self.DOWNLOAD_COMPLETE]]['file accession'])
        # compute the number of files where the specified column is True and that are valid to download
        else:
            return len(filtered[filtered[self.DOWNLOAD_VALID]]['file accession'])

    def get_file_list(self, key: str) -> set[str]:
        """
        Filters the metadata down to entries where the specified download column is True,
        and then return a sorted list of associated file paths.
        :param key: name of one of the boolean download columns
        :return: sorted list of filtered file paths
        """
        return set([x[0] for x in self._get_filtered_download_frame(key)[['file path']].values])

    def get_failed_ids(self) -> set[str]:
        """
        Filters the metadata down to entries that failed to download.
        :return: IDs of files that failed to download
        """
        filtered = self._get_filtered_download_frame(self.DOWNLOAD_FAILED)['file accession']

        return set(filtered.values)

    def get_download_list(self) -> list[tuple[str, str, str, str]]:
        """
        Generates a list of valid files to download.
        :return: sorted list of file accession, download URL, download path and download directory
        """
        # filter the metadata and extract the relevant columns
        filtered = self._get_filtered_download_frame(self.DOWNLOAD_LIST)
        filtered = filtered[['file accession', 'file download url', 'file path', 'file directory']]

        return sorted([tuple(x) for x in filtered.values])

    def add_failed_downloads(self, failed_ids: set[str]):
        """
        Increments the failed counter of the specified file entries by 1.
        :param failed_ids: set of file accessions
        """
        self.df.loc[self.df['file accession'].isin(failed_ids), self.DOWNLOAD_FAILED] += 1

    def to_tsv(self, tsv_path: str):
        """
        Writes the sorted metadata in a tab-separated file with column names in the original order.
        :param tsv_path: path to the tab separated file
        """
        self.df.to_csv(tsv_path, sep='\t', index=False)

    def to_db(self, db_path: str):
        """
        Writes the metadata in a SQLite database file.
        :param db_path: path to the SQLite database file
        """
        connection = sqlite3.connect(db_path)
        self.df.to_sql('data', connection, if_exists='replace', index=False)
        connection.close()

    def filter_download(self):
        """
        Only keep entries of files that are valid and were downloaded successfully.
        """
        self.df = self.df[(self.df[Metadata.DOWNLOAD_COMPLETE]) & (self.df[Metadata.DOWNLOAD_VALID])]

    def update_files(self):
        """
        Update the metadata.tsv and metadata.db with the current information.
        """
        self.to_tsv(self.tsv_path)
        self.to_db(self.db_path)
