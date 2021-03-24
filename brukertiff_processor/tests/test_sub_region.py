#!/usr/bin/env python
import os

import pytest
from base_processor.tests import init_ssm
from base_processor.tests import setup_processor
# our module(s)
from brukertiff_processor import BRUKERTIFFProcessor
from moto import mock_s3
from moto import mock_ssm

test_processor_data = [
    'TSeries-06272017-1624-103.brukertiff.gz'
]


@pytest.mark.parametrize("filename", test_processor_data)
def test_brukertiff_parallel_processor(filename):
    mock_ssm().start()
    mock_s3().start()
    init_ssm()
    # init task
    for x in range(2):
        for y in range(2):
            for z in range(2):
                for c in range(2):
                    for t in range(2):
                        # Create sub_region file
                        open('sub_x_{}_2_y_{}_2_z_{}_2_c_{}_2_t_{}_2.txt'.format(x, y, z, c, t), 'w').write("")
                        # Setup task with inputs
                        inputs = {
                            'file': os.path.join('/test-resources', filename),
                            'sub_region_file': 'sub_x_{}_2_y_{}_2_z_{}_2_c_{}_2_t_{}_2.txt'.format(x, y, z, c, t)
                        }
                        task = BRUKERTIFFProcessor(inputs=inputs)
                        setup_processor(task)
                        # run
                        task.run()
    mock_s3().stop()
    mock_ssm().stop()
