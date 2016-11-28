import mock
import pytest

from bcbio.distributed import objectstore
from bcbio.distributed.objectstore import GoogleDrive, GoogleDownloader


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


def test_creates_google_credentials(mock_api):
    GoogleDrive()
    objectstore.ServiceAccountCredentials.from_json_keyfile_name\
        .assert_called_once_with(
            GoogleDrive.GOOGLE_API_KEY_FILE, scopes=GoogleDrive.SCOPES)


def test_api_scope_includes_google_drive(mock_api):
    drive_scope = 'https://www.googleapis.com/auth/drive'
    assert drive_scope in GoogleDrive.SCOPES


def test_filename_with_json_key_is_present(mock_api):
    assert GoogleDrive.GOOGLE_API_KEY_FILE
    assert GoogleDrive.GOOGLE_API_KEY_FILE.endswith('.json')


def test_creates_http_auth(mock_api):
    Credentials = objectstore.ServiceAccountCredentials
    GoogleDrive()
    objectstore.build.assert_called_once_with(
        'drive', 'v3',
        Credentials.from_json_keyfile_name.return_value
    )


def test_has_a_service_attribute(mock_api):
    drive = GoogleDrive()
    assert drive.service == objectstore.build.return_value


def test_can_load_file_by_id(mock_api):
    drive = GoogleDrive()
    output_file = 'test_file'
    file_id = 'test_file_id'
    drive.download_file(file_id, output_file)
    drive.service.files().get_media.assert_called_once_with(fileId=file_id)


def test_downloader(mock_api):
    downloader = GoogleDownloader()
    assert downloader


def test_downloader_executes_request(mock_api):
    downloader = GoogleDownloader()
    fd, request = mock.Mock(), mock.Mock()
    downloader.load_to_file(fd, request)
    objectstore.http.MediaIoBaseDownload.assert_called_once_with(
        fd, request, chunksize=GoogleDownloader.CHUNK_SIZE)
