from __future__ import print_function

from datetime import datetime
import functools
import os
import pytz
import re

import matplotlib
matplotlib.use('Agg')
import pylab
pylab.rcParams['figure.figsize'] = (35.0, 12.0)
import pandas as pd

from bcbio import utils
from bcbio.graph.collectl import load_collectl
# from bcbiovm.graph.elasticluster import fetch_collectl


def get_bcbio_timings(path):
    """Fetch timing information from a bcbio log file."""
    with open(path, 'r') as fp:
        steps = {}
        for line in fp:
            matches = re.search(r'^\[([^\]]+)\] ([^:]+: .*)', line)
            if not matches:
                continue

            tstamp = matches.group(1)
            msg = matches.group(2)

            if not msg.find('Timing: ') >= 0:
                continue

            when = datetime.strptime(tstamp, '%Y-%m-%dT%H:%MZ').replace(
                tzinfo=pytz.timezone('UTC'))
            step = msg.split(":")[-1].strip()
            steps[when] = step

        return steps


def this_and_prev(iterable):
    """Walk an iterable, returning the current and previous items
    as a two-tuple."""
    try:
        item = next(iterable)
        while True:
            next_item = next(iterable)
            yield item, next_item
            item = next_item
    except StopIteration:
        return


def delta_from_prev(prev_values, tstamps, value):
    try:
        prev_val = next(prev_values)
        cur_tstamp, prev_tstamp = next(tstamps)
    except StopIteration:
        return 0

    # Take the difference from the previous value and divide by the interval
    # since the previous sample, so we always return values in units/second.
    return (value - prev_val) / (cur_tstamp - prev_tstamp).seconds


def calc_deltas(df, series=[]):
    """Many of collectl's data values are cumulative (monotonically
    increasing), so subtract the previous value to determine the value
    for the current interval.
    """
    df = df.sort(ascending=False)

    for s in series:
        prev_values = iter(df[s])
        # Burn the first value, so the first row we call delta_from_prev()
        # for gets its previous value from the second row in the series,
        # and so on.
        next(prev_values)
        df[s] = df[s].apply(functools.partial(
            delta_from_prev, iter(prev_values),
            this_and_prev(iter(df.index))))

    return df


def remove_outliers(series, stddev):
    """Remove the outliers from a series."""
    return series[(series - series.mean()).abs() < stddev * series.std()]


def prep_for_graph(df, series=[], delta_series=[], smoothing=None,
                   outlier_stddev=None):
    """Prepare a dataframe for graphing by calculating deltas for
    series that need them, resampling, and removing outliers.
    """
    graph = calc_deltas(df, delta_series)

    for s in series + delta_series:
        if smoothing:
            graph[s] = graph[s].resample(smoothing)
        if outlier_stddev:
            graph[s] = remove_outliers(graph[s], outlier_stddev)

    return graph[series + delta_series]


def add_common_plot_features(plot, steps):
    """Add plot features common to all plots, such as bcbio step
    information.
    """
    plot.yaxis.set_tick_params(labelright=True)
    plot.set_xlabel('')

    ymax = plot.get_ylim()[1]
    ticks = {}
    for tstamp, step in steps.iteritems():
        if step == 'finished':
            continue
        plot.vlines(tstamp, 0, ymax, linestyles='dashed')
        ticks[tstamp] = step
    tick_kvs = sorted(ticks.iteritems())
    top_axis = plot.twiny()
    top_axis.set_xlim(*plot.get_xlim())
    top_axis.set_xticks([k for k, v in tick_kvs])
    top_axis.set_xticklabels([v for k, v in tick_kvs], rotation=45, ha='left', size=16)
    plot.set_ylim(0)

    return plot


def graph_cpu(df, steps, num_cpus):
    graph = prep_for_graph(
        df, delta_series=['cpu_user', 'cpu_sys', 'cpu_wait'])

    graph['cpu_user'] /= 100.0
    graph['cpu_sys'] /= 100.0
    graph['cpu_wait'] /= 100.0

    plot = graph.plot()
    plot.set_ylabel('CPU core usage')
    plot.set_ylim(0, num_cpus)
    add_common_plot_features(plot, steps)

    return plot


def graph_net_bytes(df, steps, ifaces):
    series = []
    for iface in ifaces:
        series.extend(['{}_rbyte'.format(iface), '{}_tbyte'.format(iface)])

    graph = prep_for_graph(df, delta_series=series)

    for iface in ifaces:
        old_series = '{}_rbyte'.format(iface)
        new_series = '{}_receive'.format(iface)
        graph[new_series] = graph[old_series] * 8 / 1024 / 1024
        del graph[old_series]

        old_series = '{}_tbyte'.format(iface)
        new_series = '{}_transmit'.format(iface)
        graph[new_series] = graph[old_series] * 8 / 1024 / 1024
        del graph[old_series]

    plot = graph.plot()
    plot.set_ylabel('mbits/s')
    add_common_plot_features(plot, steps)

    return plot


def graph_net_pkts(df, steps, ifaces):
    series = []
    for iface in ifaces:
        series.extend(['{}_rpkt'.format(iface), '{}_tpkt'.format(iface)])

    graph = prep_for_graph(df, delta_series=series)

    plot = graph.plot()
    plot.set_ylabel('packets/s')
    add_common_plot_features(plot, steps)

    return plot


def graph_memory(df, steps, total_mem):
    graph = prep_for_graph(
        df, series=['mem_total', 'mem_free', 'mem_buffers', 'mem_cached'])

    free_memory = graph['mem_free'] + graph['mem_buffers'] + \
        graph['mem_cached']
    graph = (graph['mem_total'] - free_memory) / 1024 / 1024

    plot = graph.plot()
    plot.set_ylabel('gbytes')
    plot.set_ylim(0, total_mem)
    add_common_plot_features(plot, steps)

    return plot


def graph_disk_io(df, steps, disks):
    series = []
    for disk in disks:
        series.extend([
            '{}_sectors_read'.format(disk),
            '{}_sectors_written'.format(disk),
        ])

    graph = prep_for_graph(df, delta_series=series, outlier_stddev=2)

    for disk in disks:
        old_series = '{}_sectors_read'.format(disk)
        new_series = '{}_read'.format(disk)
        graph[new_series] = graph[old_series] * 512 / 1024 / 1024
        del graph[old_series]

        old_series = '{}_sectors_written'.format(disk)
        new_series = '{}_write'.format(disk)
        graph[new_series] = graph[old_series] * 512 / 1024 / 1024
        del graph[old_series]

    plot = graph.plot()
    plot.set_ylabel('mbytes/s')
    add_common_plot_features(plot, steps)

    return plot


def generate_graphs(collectl_datadir, bcbio_log_path, outdir, verbose=False):
    """Generate all graphs for a bcbio run."""
    if verbose:
        print('Reading timings from bcbio log...')
    steps = get_bcbio_timings(bcbio_log_path)
    start_time = min(steps.keys())
    end_time = max(steps.keys())

    if verbose:
        print('Parsing performance data...')

    dfs = {}
    hardware_info = {}
    for item in sorted(os.listdir(collectl_datadir)):
        if not item.endswith('.raw.gz'):
            continue

        df, hardware = load_collectl(
            os.path.join(collectl_datadir, item), start_time, end_time)
        if len(df) == 0:
            continue

        host = re.sub(r'-\d{8}-\d{6}\.raw\.gz$', '', item)
        hardware_info[host] = hardware
        if host not in dfs:
            dfs[host] = df
        else:
            old_df = dfs[host]
            dfs[host] = pd.concat([old_df, df])

    for host, df in dfs.iteritems():
        if verbose:
            print('Generating CPU graph for {}...'.format(host))
        graph = graph_cpu(df, steps, hardware_info[host]['num_cpus'])
        graph.get_figure().savefig(
            os.path.join(outdir, '{}_cpu.png'.format(host)),
            bbox_inches='tight', pad_inches=0.25)
        pylab.close()

        ifaces = set([
            series.split('_')[0]
            for series
             in df.keys()
             if series.startswith(('eth', 'ib'))
        ])

        if verbose:
            print('Generating network graphs for {}...'.format(host))
        graph = graph_net_bytes(df, steps, ifaces)
        graph.get_figure().savefig(
            os.path.join(outdir, '{}_net_bytes.png'.format(host)),
            bbox_inches='tight', pad_inches=0.25)
        pylab.close()

        graph = graph_net_pkts(df, steps, ifaces)
        graph.get_figure().savefig(
            os.path.join(outdir, '{}_net_pkts.png'.format(host)),
            bbox_inches='tight', pad_inches=0.25)
        pylab.close()

        if verbose:
            print('Generating memory graph for {}...'.format(host))
        graph = graph_memory(df, steps, hardware_info[host]["memory"])
        graph.get_figure().savefig(
            os.path.join(outdir, '{}_memory.png'.format(host)),
            bbox_inches='tight', pad_inches=0.25)
        pylab.close()

        if verbose:
            print('Generating storage I/O graph for {}...'.format(host))
        drives = set([
            series.split('_')[0]
            for series
             in df.keys()
             if series.startswith(('sd', 'vd', 'hd', 'xvd'))
        ])
        graph = graph_disk_io(df, steps, drives)
        graph.get_figure().savefig(
            os.path.join(outdir, '{}_disk_io.png'.format(host)),
            bbox_inches='tight', pad_inches=0.25)
        pylab.close()


def add_subparser(subparsers):
    parser = subparsers.add_parser("graph",
                                   help="Generate system graphs "
                                        "(CPU/memory/network/disk I/O "
                                        "consumption) from bcbio runs")
    parser.add_argument("log",
                        help="Local path to bcbio log file written by the run.")
    parser.add_argument("-o", "--outdir", default="monitoring/graphs",
                        help="Directory to write graphs to.")
    parser.add_argument("-r", "--rawdir", default="monitoring/collectl", required=True,
                        help="Directory to put raw collectl data files.")
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="Emit verbose output")

    return parser


def bootstrap(args):
    generate_graphs(args.rawdir, args.log, utils.safe_makedir(args.outdir),
        verbose=args.verbose)
