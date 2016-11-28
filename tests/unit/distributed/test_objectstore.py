import pytest

from bcbio.distributed import objectstore
from bcbio.distributed.objectstore import GoogleDrive


@pytest.fixture
def mock_api(mocker):
    mocker.patch('bcbio.distributed.objectstore.ServiceAccountCredentials')
    mocker.patch('bcbio.distributed.objectstore.Http')
    mocker.patch('bcbio.distributed.objectstore.build')
    mocker.patch('bcbio.distributed.objectstore.http')
    yield None


def test_create_google_drive_service(mock_api):
    service = GoogleDrive()
    assert service
