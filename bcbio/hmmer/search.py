"""Search using HMMER's online REST interface.

http://hmmer.janelia.org/search

From Nick Loman's EntrezAjax:

https://github.com/nickloman/entrezajax
"""
from __future__ import print_function
from six.moves import urllib
import logging

class SmartRedirectHandler(urllib.request.HTTPRedirectHandler):
    def http_error_302(self, req, fp, code, msg, headers):
        logging.debug(headers)
        return headers

def _hmmer(endpoint, args1, args2):
    opener = urllib.request.build_opener(SmartRedirectHandler())
    urllib.request.install_opener(opener);

    params = urllib.parse.urlencode(args1)
    try:
        req = urllib.request.Request(endpoint,
                              data = params,
                              headers={"Accept" : "application/json"})
        v = urllib.request.urlopen(req)
    except urllib.error.HTTPError as e:
        raise Exception("HTTP Error 400: %s" % e.read())

    results_url = v['location']

    enc_res_params = urllib.parse.urlencode(args2)
    modified_res_url = results_url + '?' + enc_res_params

    results_request = urllib.request.Request(modified_res_url)
    f = urllib.request.urlopen(results_request)
    return f

def phmmer(**kwargs):
    """Search a protein sequence against a HMMER sequence database.

    Arguments:
      seq - The sequence to search -- a Fasta string.
      seqdb -- Sequence database to search against.
      range -- A string range of results to return (ie. 1,10 for the first ten)
      output -- The output format (defaults to JSON).
    """
    logging.debug(kwargs)
    args = {'seq' : kwargs.get('seq'),
            'seqdb' : kwargs.get('seqdb')}
    args2 = {'output' : kwargs.get('output', 'json'),
             'range' : kwargs.get('range')}
    return _hmmer("http://hmmer.janelia.org/search/phmmer", args, args2)

def hmmscan(**kwargs):
    logging.debug(kwargs)
    args = {'seq' : kwargs.get('seq'),
            'hmmdb' : kwargs.get('hmmdb')}
    args2 = {'output' : 'json'}
    range = kwargs.get('range', None)
    if range:
        args2['range'] = range
    return _hmmer("http://hmmer.janelia.org/search/hmmscan", args, args2)

def test():
    seq = """>lcl||YPD4_1219|ftsK|128205128 putative cell division protein
MSQEYTEDKEVTLKKLSNGRRLLEAVLIVVTILAAYLMVALVSFNPSDPSWSQTAWHEPI
HNLGGSIGAWMADTLFSTFGVLAYAIPPIMVIFCWTAFRQRDASEYLDYFALSLRLIGTL
ALILTSCGLAALNIDDLYYFASGGVIGSLFSNAMLPWFNGVGATLTLLCIWVVGLTLFTG
WSWLVIAEKIGAAVLGSLTFITNRSRREERYDDEDSYHDDDHADGRDITGQEKGVVSNKG
VVSNNAVVGAGVAASSALAHGDDDVLFSAPSVTDSIVEHGSVVATGTETTDTKATDTNDE
YDPLLSPLRATDYSVQDATSSPIADVAVEPVLNHDAAAIYGTTPVMTNTATPPLYSFELP
EESLPIQTHAAPTERPEPKLGAWDMSPTPVSHSPFDFSAIQRPVGQLESRQPGSNQSGSH
QIHSAQSSHISVGNTPYMNPGLDAQIDGLSTTSLTNKPVLASGTVAAATAAAAFMPAFTA
TSDSSSQIKQGIGPELPRPNPVRIPTRRELASFGIKLPSQRMAEQELRERDGDETQNPQM
AASSYGTEITSDEDAALQQAILRKAFADQQSERYALSTLAEQSSITERSPAAEMPTTPSQ
VSDLEDEQALQEAELRQAFAAQQQHRYGATGDTDNAVDNIRSVDTSTAFTFSPIADLVDD
SPREPLFTLSPYVDETDVDEPVQLEGKEESLLQDYPEQVPTYQPPVQQAHLGQSAPTQPS
HTQSTYGQSTYGQSTYGQSTPAPVSQPVVTSASAISTSVTPTSIASLNTAPVSAAPVAPS
PQPPAFSQPTAAMDSLIHPFLMRNDQPLQKPTTPLPTLDLLSSPPAEEEPVDMFALEQTA
RLVEARLGDYRVKAEVVGISPGPVITRFELDLAPGVKASRISNLSRDLARSLSAIAVRVV
EVIPGKPYVGLELPNKHRQTVYLREVLDCAKFRENPSPLAIVLGKDIAGQPVVADLAKMP
HLLVAGTTGSGKSVGVNAMILSILYKATPDDVRFIMIDPKMLELSVYEGIPHLLTGVVTD
MKDAANALRWCVGEMERRYKLMSALGVRNLAGYNERVAQAEAMGRPIPDPFWKPSDSMDI
SPPMLVKLPYIVVMVDEFADLMMTVGKKVEELIARLAQKARAAGIHLVLATQRPSVDVIT
GLIKANIPTRIAFTVSSKIDSRTILDQGGAESLLGMGDMLYMAPNSSIPVRVHGAFVRDQ
EVHAVVNDWKARGRPQYIDSILSGGEEGEGGGLGLDSDEELDPLFDQAVNFVLEKRRASI
SGVQRQFRIGYNRAARIIEQMEAQQIVSTPGHNGNREVLAPPPHE"""
    handle = hmmscan(hmmdb = 'pfam', seq = seq)
    import json
    j = json.loads(handle.read())
    print(json.dumps(j, sort_keys=True, indent=4))

# test()

