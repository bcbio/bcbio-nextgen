import os
import unittest
from types import StringType

from bcbio import broad
from bcbio.pipeline import config_utils

from nose.plugins.attrib import attr

'''
To call this test individually, 
python test_picard_version.py ValidatePicard.test_version
'''

class ValidatePicard(unittest.TestCase):
    def setUp(self):
        self._config_file = os.path.join(os.path.dirname(__file__), 
                                         '../../config', 
                                         'bcbio_system.yaml' )
        config, config_file = config_utils.load_system_config(self._config_file)
        self._picard = broad.runner_from_config(config)

    @attr("unit")    
    def test_version(self):
        '''
        Before we call picard, we check the verison. Do we parse it
        correctly?  see "Version checking for Picard ViewSam not
        robust to environmental java variables #494" for more detail.
        '''
        version = self._picard.get_picard_version("ViewSam")
        assert type(version) is StringType, "version is not a string: {0}".format(version)


if __name__ == "__main__":
    unittest.main()
