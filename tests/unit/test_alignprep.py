from bcbio.ngsalign import alignprep


class TestMisc(object):
    """Additional unit test cases to run regularly to confirm code logic.
    """
    def test_align_split_size(self):
        """Checks on logic for estimating align split size.
        """
        assert alignprep._pick_align_split_size(10, 5, 20, 50) == 20
        assert alignprep._pick_align_split_size(250, 5, 20, 50) == 20
        assert alignprep._pick_align_split_size(500, 5, 20, 50) == 40
        assert alignprep._pick_align_split_size(750, 5, 20, 50) == 60
