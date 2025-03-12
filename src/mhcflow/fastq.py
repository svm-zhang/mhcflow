import subprocess as sp

import numpy as np
import polars as pl


def dump_to_fastq(qnames: pl.Series, nproc: int = 4):
    print(qnames.shape)
    qname_batches = np.array_split(qnames.to_numpy(), nproc)
    print(len(qname_batches))
