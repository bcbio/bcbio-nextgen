import calendar
import glob
import gzip
import math
import os.path
import re

import pandas as pd

def _parse_raw(fp, start_tstamp, end_tstamp):
    import progressbar
    widgets = [
        os.path.basename(fp.name), ': ',
        progressbar.Bar(marker='-', left='[', right=']'), ' ',
        progressbar.Percentage(), ' ', progressbar.ETA(),
    ]
    # We don't know what the file's uncompressed size will wind up being,
    # so take an educated guess and ignore the AssertionError later on
    # if it winds up being bigger than we guess.
    bar = progressbar.ProgressBar(
        widgets=widgets, maxval=os.path.getsize(fp.name) * 15)
    bar.start()
    bar.update(0)

    tstamp = 0
    hardware = {}
    data = {}
    for line in fp:
        matches = re.search(r'^>>> (\d+).\d+ <<<', line)
        if matches:
            try:
                bar.update(fp.tell())
            except AssertionError:
                pass

            tstamp = int(matches.group(1))
            if (tstamp >= start_tstamp) or (tstamp <= end_tstamp):
                data[tstamp] = {
                    'disk': {},
                    'mem': {},
                    'net': {},
                    'proc': {},
                }
            continue

        if line.startswith('# SubSys: '):
            matches = re.search(r'\sNumCPUs: (\d+)\s+', line)
            if matches:
                hardware['num_cpus'] = int(matches.group(1))
            continue
        if line.startswith('# Kernel: '):
            matches = re.search(r'\sMemory: (\d+)\s+kB', line)
            if matches:
                hardware['memory'] = int(math.ceil(float(matches.group(1)) / math.pow(1024.0, 2.0)))
            continue

        if (tstamp < start_tstamp) or (tstamp > end_tstamp):
            continue

        if line.startswith('cpu '):
            # Don't know what the last two fields are, but they
            # always seem to be 0, and collectl doesn't parse them
            # in formatit::dataAnalyze().
            (title, user, nice, sys, idle, wait, irq,
             soft, steal) = line.split()[:9]
            data[tstamp]['cpu'] = {
                 'user': user,
                 'nice': nice,
                 'sys': sys,
                 'idle': idle,
                 'wait': wait,
                 'irq': irq,
                 'soft': soft,
                 'steal': steal,
            }
        elif line.startswith('disk '):
            (title, major, minor, node,
             num_reads, reads_merged, sectors_read, msec_spent_reading,
             num_writes, writes_merged, sectors_written, msec_spent_writing,
             iops_in_progress, msec_spent_on_iops,
             weighted_msec_spent_on_iops) = line.split()
            data[tstamp]['disk'][node] = {
                'num_reads': num_reads,
                'reads_merged': reads_merged,
                'sectors_read': sectors_read,
                'msec_spent_reading': msec_spent_reading,
                'num_writes': num_writes,
                'writes_merged': writes_merged,
                'sectors_written': sectors_written,
                'msec_spent_writing': msec_spent_writing,
                'iops_in_progress': iops_in_progress,
                'msec_spent_on_iops': msec_spent_on_iops,
                'weighted_msec_spent_on_iops': weighted_msec_spent_on_iops,
            }
        elif line.startswith('Net '):
            # Older kernel versions don't have whitespace after
            # the interface colon:
            #
            #   Net   eth0:70627391
            #
            # unlike newer kernels:
            #
            #   Net   eth0: 415699541
            line = re.sub(r'^(Net\s+[^:]+):', r'\1: ', line)

            (title, iface,
             rbyte, rpkt, rerr, rdrop, rfifo,
             rframe, rcomp, rmulti,
             tbyte, tpkt, terr, tdrop, tfifo,
             tcoll, tcarrier, tcomp) = line.split()
            iface = iface.replace(':', '')
            data[tstamp]['net'][iface] = {
                 'rbyte': rbyte,
                 'rpkt': rpkt,
                 'rerr': rerr,
                 'rdrop': rdrop,
                 'rfifo': rfifo,
                 'rframe': rframe,
                 'rcomp': rcomp,
                 'rmulti': rmulti,
                 'tbyte': tbyte,
                 'tpkt': tpkt,
                 'terr': terr,
                 'tdrop': tdrop,
                 'tfifo': tfifo,
                 'tcoll': tcoll,
                 'tcarrier': tcarrier,
                 'tcomp': tcomp,
            }
        elif line.startswith('MemTotal:'):
            title, amount, unit = line.split()
            data[tstamp]['mem']['total'] = amount
        elif line.startswith('MemFree:'):
            title, amount, unit = line.split()
            data[tstamp]['mem']['free'] = amount
        elif line.startswith('Buffers:'):
            title, amount, unit = line.split()
            data[tstamp]['mem']['buffers'] = amount
        elif line.startswith('Cached:'):
            title, amount, unit = line.split()
            data[tstamp]['mem']['cached'] = amount
        # We don't currently do anything with process data,
        # so don't bother parsing it.
        elif False and line.startswith('proc:'):
            title_pid, rest = line.split(None, 1)
            title, pid = title_pid.split(':')

            if pid not in data[tstamp]['proc']:
                data[tstamp]['proc'][pid] = {}

            if rest.startswith('cmd '):
                title, cmd = rest.split(None, 1)
                data[tstamp]['proc'][pid]['cmd'] = cmd
            elif rest.startswith('io read_bytes: '):
                value = rest.split(':')[1].strip()
                data[tstamp]['proc'][pid]['read_bytes'] = value
            elif rest.startswith('io write_bytes: '):
                value = rest.split(':')[1].strip()
                data[tstamp]['proc'][pid]['write_bytes'] = value

    bar.finish()
    return hardware, data


class _CollectlGunzip(gzip.GzipFile):
    """collectl writes data to its files incrementally, and doesn't
    add a CRC to the end until it rotates the log. Ignore the CRC
    errors; they're innocuous in this case.
    """

    def _read_eof(self):
        return


def load_collectl(pattern, start_time, end_time):
    """Read data from collectl data files into a pandas DataFrame.
        :pattern: Absolute path to raw collectl files
    """
    start_tstamp = calendar.timegm(start_time.utctimetuple())
    end_tstamp = calendar.timegm(end_time.utctimetuple())

    cols = []
    rows = []

    for path in glob.glob(pattern):
        hardware, raw = _parse_raw(
            _CollectlGunzip(path, 'r'), start_tstamp, end_tstamp)

        if not cols:
            instances = {
                'disk': set(),
                'net': set(),
                'proc': set(),
            }
            for tstamp, sample in raw.items():
                for group, items in sample.items():
                    if group == 'disk':
                        instances['disk'] = instances['disk'].union(
                            items.keys())
                    elif group == 'net':
                        instances['net'] = instances['net'].union(
                            items.keys())
                    elif group == 'proc':
                        instances['proc'] = instances['proc'].union(
                            items.keys())

            cols = ['tstamp']
            cols.extend([
                'cpu_{}'.format(var)
                for var
                 in ['user', 'nice', 'sys', 'idle', 'wait',
                     'irq', 'soft', 'steal']
                 ])
            for node in instances['disk']:
                cols.extend([
                    '{}_{}'.format(node, var)
                    for var
                     in ['num_reads', 'reads_merged',
                         'sectors_read', 'msec_spent_reading',
                         'num_writes', 'writes_merged',
                         'sectors_written', 'msec_spent_writing',
                         'iops_in_progress', 'msec_spent_on_iops',
                         'weighted_msec_spent_on_iops']
                     ])
            cols.extend([
                'mem_{}'.format(var)
                for var
                 in ['total', 'free', 'buffers', 'cached']
                 ])
            for iface in instances['net']:
                cols.extend([
                    '{}_{}'.format(iface, var)
                    for var
                     in ['rbyte', 'rpkt', 'rerr', 'rdrop',
                         'rfifo', 'rframe', 'rcomp', 'rmulti',
                         'tbyte', 'tpkt', 'terr', 'tdrop',
                         'tfifo', 'tcoll', 'tcarrier', 'tcomp']
                     ])
            for pid in instances['proc']:
                cols.extend([
                    '{}_{}'.format(pid, var)
                    for var
                     in ['name', 'read_bytes', 'write_bytes']
                     ])

        for tstamp, sample in raw.items():
            if ('cpu' not in sample or
                'disk' not in sample or
                'mem' not in sample):
                # Skip incomplete samples; there might be a truncated
                # sample on the end of the file.
                continue

            values = [tstamp]
            values.extend([
                sample['cpu']['user'], sample['cpu']['nice'],
                sample['cpu']['sys'], sample['cpu']['idle'],
                sample['cpu']['wait'], sample['cpu']['irq'],
                sample['cpu']['soft'], sample['cpu']['steal'],
            ])
            for node in instances['disk']:
                data = sample['disk'].get(node, {})
                values.extend([
                    data.get('num_reads', 0),
                    data.get('reads_merged', 0),
                    data.get('sectors_read', 0),
                    data.get('msec_spent_reading', 0),
                    data.get('num_writes', 0),
                    data.get('writes_merged', 0),
                    data.get('sectors_written', 0),
                    data.get('msec_spent_writing', 0),
                    data.get('iops_in_progress', 0),
                    data.get('msec_spent_on_iops', 0),
                    data.get('weighted_msec_spent_on_iops', 0),
                ])
            values.extend([
                sample['mem']['total'], sample['mem']['free'],
                sample['mem']['buffers'], sample['mem']['cached'],
            ])
            for iface in instances['net']:
                data = sample['net'].get(iface, {})
                values.extend([
                    data.get('rbyte', 0), data.get('rpkt', 0),
                    data.get('rerr', 0), data.get('rdrop', 0),
                    data.get('rfifo', 0), data.get('rframe', 0),
                    data.get('rcomp', 0), data.get('rmulti', 0),
                    data.get('tbyte', 0), data.get('tpkt', 0),
                    data.get('terr', 0), data.get('tdrop', 0),
                    data.get('tfifo', 0), data.get('tcoll', 0),
                    data.get('tcarrier', 0), data.get('tcomp', 0),
                ])
            if 'proc' in sample:
                for pid in instances['proc']:
                    data = sample['proc'].get(pid, {})
                    values.extend([
                        data.get('cmd', ''),
                        data.get('read_bytes', 0),
                        data.get('write_bytes', 0),
                    ])

            rows.append(values)

    if len(rows) == 0:
        return pd.DataFrame(columns=cols), {}

    df = pd.DataFrame(rows, columns=cols)
    df = df.convert_objects(convert_numeric=True)
    df['tstamp'] = df['tstamp'].astype('datetime64[s]')
    df.set_index('tstamp', inplace=True)
    df = df.tz_localize('UTC')

    return df, hardware
