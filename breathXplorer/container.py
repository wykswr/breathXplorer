from pathlib import Path
from typing import Iterator, Union

import numpy as np
import pandas as pd


class Container:
    data: pd.DataFrame

    def __init__(self, data):
        self.data = data

    @property
    def table(self):
        return self.data

    @property
    def mz(self):
        return self.data.index.values.astype(float)

    def __len__(self):
        return self.data.shape[0]

    def to_csv(self, file: Union[str, Path]):
        self.data.to_csv(file)


class FeatureSet(Container):
    @property
    def scan_time(self):
        """
        Get the scan time of the feature set.
        :return: scan time
        """
        return self.data.columns[1:].values.astype(float)

    def __getitem__(self, key):
        return self.data.loc[key, 1:].values.astype(float)

    @property
    def intensity(self):
        return self.data['intensity'].values.astype(float)


class Sample(Container):
    @property
    def sample_name(self):
        return self.data.columns

    def __getitem__(self, key):
        return self.data.loc[key, :].values.astype(float)


class TandemMS:
    feature: np.ndarray
    spectra: dict

    def __init__(self, feature):
        self.feature = feature

    def build(self, mzML: Iterator, radium: float):
        """
        Build the tandem MS spectra using the feature list.
        :param mzML: Iterator of mzML file
        :param radium: tolerance of the mz
        :return: None
        """
        self.spectra = dict()
        if self.feature.size == 0:
            raise ValueError('Empty feature')
        for sp in mzML:
            if sp['ms level'] != 2:
                continue
            pm = sp['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
            indicator = np.abs(self.feature - pm) < radium
            if indicator.sum() == 0:
                continue
            self.spectra[pm] = sp
        self.spectra = dict(sorted(self.spectra.items(), key=lambda x: x[0]))

    @staticmethod
    def __translate(mz: np.ndarray, intensity: np.ndarray, precursor: float) -> str:
        return 'BEGIN IONS\n' + \
            f'PEPMASS={precursor}\n' + \
            'MSLEVEL=2\n' + \
            'CHARGE=1+\n' + \
            "\n".join([f"{mz[i]} {intensity[i]}" for i in range(mz.size) if intensity[i] > 0.001]) + \
            '\nEND IONS\n\n'

    def to_mgf(self, file: Union[str, Path]):
        """
        Write the spectra to an MGF file.
        :param file: MGF file path
        :return: None
        """
        with open(file, 'w') as handle:
            for pm, sp in self.spectra.items():
                mz, intensity = sp['m/z array'], sp['intensity array']
                handle.write(self.__translate(mz, intensity, pm))
