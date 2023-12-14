from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterator, Union

import numpy as np
import pandas as pd
from pyteomics import mzml, mzxml

from .utils import annotate_adduct, annotate_isotope


class Spectra(ABC):
    def __init__(self, sp):
        self.sp = sp

    @property
    @abstractmethod
    def mz(self) -> np.ndarray:
        pass

    @property
    @abstractmethod
    def intensity(self) -> np.ndarray:
        pass

    @property
    @abstractmethod
    def level(self) -> int:
        pass

    @property
    @abstractmethod
    def tic(self) -> float:
        pass

    @property
    @abstractmethod
    def scan_start_time(self) -> float:
        pass

    @property
    @abstractmethod
    def precursor(self) -> float:
        pass


class MzmlSpectra(Spectra):
    @property
    def mz(self):
        return self.sp['m/z array']

    @property
    def intensity(self):
        return self.sp['intensity array']

    @property
    def level(self):
        return self.sp['ms level']

    @property
    def tic(self):
        return self.sp['total ion current']

    @property
    def scan_start_time(self):
        return self.sp['scanList']['scan'][0]['scan start time']

    @property
    def precursor(self):
        return self.sp['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']


class MzxmlSpectra(Spectra):
    @property
    def mz(self):
        return self.sp["m/z array"]

    @property
    def intensity(self):
        return self.sp["intensity array"]

    @property
    def level(self):
        return self.sp["msLevel"]

    @property
    def tic(self):
        return self.sp["totIonCurrent"]

    @property
    def scan_start_time(self):
        return self.sp["retentionTime"]

    @property
    def precursor(self):
        return self.sp["precursorMz"][0]["precursorMz"]


class Source:
    def __init__(self, file: Union[str, Path]):
        self.file = Path(file)

    def __iter__(self) -> Iterator[Spectra]:
        if self.file.suffix == '.mzML':
            with mzml.read(str(self.file.absolute())) as handle:
                for sp in handle:
                    yield MzmlSpectra(sp)
        elif self.file.suffix == '.mzXML':
            with mzxml.read(str(self.file.absolute())) as handle:
                for sp in handle:
                    yield MzxmlSpectra(sp)
        else:
            raise ValueError('Unsupported file format')


class Container:
    data: pd.DataFrame

    def __init__(self, data):
        self.data = data

    @property
    def table(self) -> pd.DataFrame:
        return self.data

    @property
    def mz(self) -> np.ndarray:
        return self.data.index.values.astype(float)

    def __len__(self):
        return self.data.shape[0]

    def to_csv(self, file: Union[str, Path], adduct: bool = False, isotope: bool = False):
        # deep copy
        data = self.data.copy()
        data = data.applymap(lambda x: np.round(x, 0))
        features = data.index.values.astype(float)
        features = np.round(features, 4)
        data.insert(0, 'm/z', features)
        data.insert(0, 'ID', np.arange(data.shape[0]))
        if adduct:
            data.insert(1, 'adduct', annotate_adduct(features, .001))
        if isotope:
            data.insert(1, 'isotope', annotate_isotope(features, .001))
        data.to_csv(file, index=False)


class FeatureSet(Container):
    @property
    def scan_time(self) -> np.ndarray:
        """
        Get the scan time of the feature set.
        :return: scan time
        """
        return self.data.columns[1:].values.astype(float)

    def __getitem__(self, key: float) -> np.ndarray:
        return self.data.loc[key, 1:].values.astype(float)

    @property
    def rsd(self) -> pd.Series:
        return self.data.iloc[:, 1:].apply(lambda x: np.std(x) / np.mean(x), axis=1)

    @property
    def intensity(self) -> np.ndarray:
        return self.data['intensity'].values.astype(float)

    def rsd_control(self, threshold: float) -> 'FeatureSet':
        """
        Filter the feature set by RSD.
        :param threshold: threshold of RSD
        :return: filtered feature set
        """
        return FeatureSet(self.data[self.rsd > threshold])


class Sample(Container):
    @property
    def sample_name(self):
        return self.data.columns

    def __getitem__(self, key: float):
        return self.data.loc[key, :].values.astype(float)


class TandemMS:
    feature: np.ndarray
    spectra: dict

    def __init__(self, feature):
        self.feature = feature

    def build(self, source: Source, radium: float):
        """
        Build the tandem MS spectra using the feature list.
        :param source: source of the MS data
        :param radium: tolerance of the mz
        :return: None
        """
        self.spectra = dict()
        if self.feature.size == 0:
            raise ValueError('Empty feature')
        for sp in source:
            if sp.level != 2:
                continue
            indicator = np.abs(self.feature - sp.precursor) < radium
            if indicator.sum() == 0:
                continue
            self.spectra[sp.precursor] = sp
        self.spectra = dict(sorted(self.spectra.items(), key=lambda x: x[0]))

    @staticmethod
    def __translate(mz: np.ndarray, intensity: np.ndarray, precursor: float) -> str:
        return 'BEGIN IONS\n' + \
            f'PEPMASS={precursor}\n' + \
            'MSLEVEL=2\n' + \
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
                handle.write(self.__translate(sp.mz, sp.intensity, pm))
