import mock
import pytest

from bcbio.distributed import objectstore
from bcbio.distributed.objectstore import GoogleDriveServiceFactory
from bcbio.distributed.objectstore import GoogleDownloader
from bcbio.distributed.objectstore import GoogleDrive


@pytest.fixture
def mock_api(mocker):
    mocker.patch('bcbio.distributed.objectstore.ServiceAccountCredentials')
    mocker.patch('bcbio.distributed.objectstore.Http')
    mocker.patch('bcbio.distributed.objectstore.build')
    mock_http = mocker.patch('bcbio.distributed.objectstore.http')
    mocker.patch('bcbio.distributed.objectstore.open')
    media = mock_http.MediaIoBaseDownload.return_value
    media.next_chunk.side_effect = [
        (mock.Mock(), True)
    ]
    yield None


def test_create_google_drive_service(mock_api):
    service = GoogleDriveServiceFactory.create('test')
    assert service


def test_creates_google_credentials(mock_api):
    GoogleDriveServiceFactory.create('TEST_API_KEY_FILE')
    objectstore.ServiceAccountCredentials.from_json_keyfile_name\
        .assert_called_once_with(
            'TEST_API_KEY_FILE',
            scopes=GoogleDriveServiceFactory.SCOPES
        )


def test_api_scope_includes_google_drive(mock_api):
    drive_scope = 'https://www.googleapis.com/auth/drive'
    assert drive_scope in GoogleDriveServiceFactory.SCOPES


def test_creates_api_credentials(mock_api):
    cred = objectstore.ServiceAccountCredentials.from_json_keyfile_name()
    GoogleDriveServiceFactory.create('TEST')
    objectstore.build.assert_called_once_with(
        'drive', 'v3',
        cred.authorize.return_value
    )


def test_creates_http_auth(mock_api):
    cred = objectstore.ServiceAccountCredentials.from_json_keyfile_name()
    GoogleDriveServiceFactory.create('TEST')
    cred.authorize.assert_called_once_with(objectstore.Http())


def test_downloader(mock_api):
    downloader = GoogleDownloader()
    assert downloader


def test_downloader_executes_request(mock_api):
    fd, request = mock.Mock(), mock.Mock()
    media = objectstore.http.MediaIoBaseDownload.return_value
    media.next_chunk.side_effect = [
        (mock.Mock(), True)
    ]
    downloader = GoogleDownloader()
    downloader.load_to_file(fd, request)
    objectstore.http.MediaIoBaseDownload.assert_called_once_with(
        fd, request, chunksize=GoogleDownloader.CHUNK_SIZE)


def test_loads_content_in_chunks(mock_api):
    fd, request = mock.Mock(), mock.Mock()
    media = objectstore.http.MediaIoBaseDownload.return_value
    media.next_chunk.side_effect = [
        (mock.Mock(), True)
    ]
    downloader = GoogleDownloader()
    downloader.load_to_file(fd, request)
    media.next_chunk.assert_called_once_with(
        num_retries=GoogleDownloader.NUM_RETRIES)


def test_loads_chunks_until_done(mock_api):
    fd, request = mock.Mock(), mock.Mock()
    media = objectstore.http.MediaIoBaseDownload.return_value
    next_chunk = [
        (mock.Mock(), False),
        (mock.Mock(), False),
        (mock.Mock(), True),

    ]
    media.next_chunk.side_effect = next_chunk
    downloader = GoogleDownloader()
    downloader.load_to_file(fd, request)
    assert objectstore.http.MediaIoBaseDownload().next_chunk.call_count == 3


class TestGoogleDrive(object):

    @pytest.yield_fixture
    def drive(self, mock_api):
        yield GoogleDrive()

    @pytest.mark.parametrize(('url', 'expected'), [
        ('foo.com', False),
        ('http://example.com', False),
        ('https://example.pl', False),
        ('https://drive.google.com', False),
        ('https://drive.google.com/file/d/1234', True),
        ('https://drive.google.com/file/d/1234/view', True),
    ])
    def test_check_repource(self, drive, url, expected):
        result = drive.check_resource(url)
        assert result == expected

    @pytest.mark.parametrize('url', [
        'https://drive.google.com/file/d/TEST_ID/view',
        'https://drive.google.com/file/d/TEST_ID/',
        'https://drive.google.com/file/d/TEST_ID',
    ])
    def test_parse_remote(self, drive, url):
        expected = 'TEST_ID'
        result = drive.parse_remote(url)
        assert result == expected

    def test_filename_with_json_key_is_present(self, mock_api):
        assert GoogleDrive.GOOGLE_API_KEY_FILE
        assert GoogleDrive.GOOGLE_API_KEY_FILE.endswith('.json')

    def test_can_load_file_by_id(self, drive):
        output_file = 'test_file'
        file_id = 'test_file_id'
        drive._download_file(file_id, output_file)
        drive.service.files().get_media.assert_called_once_with(fileId=file_id)

    def test_opens_output_file_for_writing(self, drive):
        drive._download_file('test_file_id', 'test_fname')
        objectstore.open.assert_called_once_with('test_fname', 'w')

    def test_downloads_file(self, drive, mocker):
        mock_load = mocker.patch.object(GoogleDownloader, 'load_to_file')
        drive._download_file('test_file_id', 'test_fname')
        fd = objectstore.open().__enter__()
        mock_load.assert_called_once_with(
            fd, drive.service.files().get_media.return_value)


