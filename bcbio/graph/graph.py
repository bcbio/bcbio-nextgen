from __future__ import print_function

from datetime import datetime
import collections
import functools
import os
import gzip
import pytz
import re
import socket

import pandas as pd
import cPickle as pickle

from bcbio import utils
from bcbio.graph.collectl import load_collectl

mpl = utils.LazyImport("matplotlib")
plt = utils.LazyImport("matplotlib.pyplot")
pylab = utils.LazyImport("pylab")

def _setup_matplotlib():
    # plt.style.use('ggplot')
    mpl.use('Agg')
    pylab.rcParams['image.cmap'] = 'viridis'
    pylab.rcParams['figure.figsize'] = (35.0, 12.0)
    # pylab.rcParams['figure.figsize'] = (100, 100)
    pylab.rcParams['figure.dpi'] = 300
    pylab.rcParams['font.size'] = 25

def get_bcbio_nodes(path):
    """Fetch the local nodes (-c local) that contain collectl files from
       the bcbio log file.

       :returns: A list with unique (non-FQDN) local hostnames
                 where collectl raw logs can be found.
    """
    with open(path, 'r') as file_handle:
        hosts = collections.defaultdict(dict)
        for line in file_handle:
            matches = re.search(r'\]\s([^:]+):', line)
            if not matches:
                continue
            # Format of the record will be "[Date] host: Timing: Step" if distributed,
            # otherwise the host will be missing and it means its a local run, we can stop
            elif 'Timing: ' in line and line.split(': ')[1] != 'Timing':
                hosts = collections.defaultdict(dict, {socket.gethostname() : {}})
                break

            hosts[matches.group(1)]

    return hosts

def get_bcbio_timings(path):
    """Fetch timing information from a bcbio log file."""
    with open(path, 'r') as file_handle:
        steps = {}
        for line in file_handle:
            matches = re.search(r'^\[([^\]]+)\] ([^:]+: .*)', line)
            if not matches:
                continue

            tstamp = matches.group(1)
            msg = matches.group(2)

            # XXX: new special logs do not have this
            #if not msg.find('Timing: ') >= 0:
            #    continue

            when = datetime.strptime(tstamp, '%Y-%m-%dT%H:%MZ').replace(
                tzinfo=pytz.timezone('UTC'))

            step = msg.split(":")[-1].strip()
            steps[when] = step

        return steps

def plot_inline_jupyter(plot):
    """ Plots inside the output cell of a jupyter notebook if %matplotlib magic
        is defined.
    """
    _setup_matplotlib()
    try:
        get_ipython()
        plt.show(plot)
    except NameError:
        pass

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
    return (prev_val - value) / (prev_tstamp - cur_tstamp).seconds


def calc_deltas(data_frame, series=None):
    """Many of collectl's data values are cumulative (monotonically
    increasing), so subtract the previous value to determine the value
    for the current interval.
    """
    series = series or []
    data_frame = data_frame.sort_index(ascending=True)

    for s in series:
        prev_values = iter(data_frame[s])
        # Burn the first value, so the first row we call delta_from_prev()
        # for gets its previous value from the second row in the series,
        # and so on.
        next(prev_values)
        data_frame[s] = data_frame[s].apply(functools.partial(
            delta_from_prev, iter(prev_values),
            this_and_prev(iter(data_frame.index))))

    return data_frame


def remove_outliers(series, stddev):
    """Remove the outliers from a series."""
    return series[(series - series.mean()).abs() < stddev * series.std()]


def prep_for_graph(data_frame, series=None, delta_series=None, smoothing=None,
                   outlier_stddev=None):
    """Prepare a dataframe for graphing by calculating deltas for
    series that need them, resampling, and removing outliers.
    """
    series = series or []
    delta_series = delta_series or []
    graph = calc_deltas(data_frame, delta_series)

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
    _setup_matplotlib()
    plot.yaxis.set_tick_params(labelright=True)
    plot.set_xlabel('')

    ymax = plot.get_ylim()[1]
    ticks = {}
    for tstamp, step in steps.iteritems():
        if step == 'finished':
            continue
        plot.vlines(tstamp, 0, ymax, linestyles='dashed')
        tstamp = mpl.dates.num2epoch(mpl.dates.date2num(tstamp))
        ticks[tstamp] = step
    tick_kvs = sorted(ticks.iteritems())
    top_axis = plot.twiny()
    top_axis.set_xlim(*plot.get_xlim())
    top_axis.set_xticks([k for k, v in tick_kvs])
    top_axis.set_xticklabels([v for k, v in tick_kvs],
                             rotation=45, ha='left', size=pylab.rcParams['font.size'])

    plot.set_ylim(0)

    return plot


def graph_cpu(data_frame, steps, num_cpus):
    graph = prep_for_graph(
        data_frame, delta_series=['cpu_user', 'cpu_sys', 'cpu_wait'])

    graph['cpu_user'] /= 100.0
    graph['cpu_sys'] /= 100.0
    graph['cpu_wait'] /= 100.0

    plot = graph.plot()
    plot.set_ylabel('CPU core usage')
    plot.set_ylim(0, num_cpus)
    add_common_plot_features(plot, steps)

    plot_inline_jupyter(plot)

    return plot, graph


def graph_net_bytes(data_frame, steps, ifaces):
    series = []
    for iface in ifaces:
        series.extend(['{}_rbyte'.format(iface), '{}_tbyte'.format(iface)])

    graph = prep_for_graph(data_frame, delta_series=series)

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
    plot.set_ylim(0, 2000)

    add_common_plot_features(plot, steps)

    plot_inline_jupyter(plot)

    return plot, graph


def graph_net_pkts(data_frame, steps, ifaces):
    series = []
    for iface in ifaces:
        series.extend(['{}_rpkt'.format(iface), '{}_tpkt'.format(iface)])

    graph = prep_for_graph(data_frame, delta_series=series)

    plot = graph.plot()
    plot.set_ylabel('packets/s')
    add_common_plot_features(plot, steps)

    plot_inline_jupyter(plot)

    return plot, graph


def graph_memory(data_frame, steps, total_mem):
    graph = prep_for_graph(
        data_frame, series=['mem_total', 'mem_free', 'mem_buffers',
                            'mem_cached'])

    free_memory = graph['mem_free'] + graph['mem_buffers'] + \
        graph['mem_cached']
    graph = (graph['mem_total'] - free_memory) / 1024 / 1024

    plot = graph.plot()
    plot.set_ylabel('gbytes')
    plot.set_ylim(0, total_mem)
    add_common_plot_features(plot, steps)

    plot_inline_jupyter(plot)

    return plot, graph


def graph_disk_io(data_frame, steps, disks):
    series = []
    for disk in disks:
        series.extend([
            '{}_sectors_read'.format(disk),
            '{}_sectors_written'.format(disk),
        ])

    graph = prep_for_graph(data_frame, delta_series=series, outlier_stddev=2)

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

    plot_inline_jupyter(plot)

    return plot, graph


def log_time_frame(bcbio_log):
    """The bcbio running time frame.

    :return:    an instance of :class collections.namedtuple:
                with the following fields: start and end
    """
    output = collections.namedtuple("Time", ["start", "end", "steps"])
    bcbio_timings = get_bcbio_timings(bcbio_log)
    return output(min(bcbio_timings), max(bcbio_timings), bcbio_timings)

def rawfile_within_timeframe(rawfile, timeframe):
    """ Checks whether the given raw filename timestamp falls within [start, end] timeframe.
    """
    matches = re.search(r'-(\d{8})-', rawfile)
    if matches:
        ftime = datetime.strptime(matches.group(1), "%Y%m%d")
        ftime = pytz.utc.localize(ftime)

    return ftime.date() >= timeframe[0].date() and ftime.date() <= timeframe[1].date()


def resource_usage(bcbio_log, cluster, rawdir, verbose):
    """Generate system statistics from bcbio runs.

    Parse the obtained files and put the information in
    a :class pandas.DataFrame:.

    :param bcbio_log:   local path to bcbio log file written by the run
    :param cluster:
    :param rawdir:      directory to put raw data files
    :param verbose:     increase verbosity

    :return: a tuple with three dictionaries, the first one contains
             an instance of :pandas.DataFrame: for each host, the second one
             contains information regarding the hardware configuration and
             the last one contains information regarding timing.
    :type return: tuple
    """
    data_frames = {}
    hardware_info = {}
    time_frame = log_time_frame(bcbio_log)

    for collectl_file in sorted(os.listdir(rawdir)):
        if not collectl_file.endswith('.raw.gz'):
            continue

        # Only load filenames within sampling timerange (gathered from bcbio_log time_frame)
        if rawfile_within_timeframe(collectl_file, time_frame):

            collectl_path = os.path.join(rawdir, collectl_file)
            data, hardware = load_collectl(
                collectl_path, time_frame.start, time_frame.end)

            if len(data) == 0:
                #raise ValueError("No data present in collectl file %s, mismatch in timestamps between raw collectl and log file?", collectl_path)
                continue

            host = re.sub(r'-\d{8}-\d{6}\.raw\.gz$', '', collectl_file)
            hardware_info[host] = hardware
            if host not in data_frames:
                data_frames[host] = data
            else:
                data_frames[host] = pd.concat([data_frames[host], data])

    return (data_frames, hardware_info, time_frame.steps)


def generate_graphs(data_frames, hardware_info, steps, outdir,
                    verbose=False):
    """Generate all graphs for a bcbio run."""
    _setup_matplotlib()
    # Hash of hosts containing (data, hardware, steps) tuple
    collectl_info = collections.defaultdict(dict)

    for host, data_frame in data_frames.iteritems():
        if verbose:
            print('Generating CPU graph for {}...'.format(host))
        graph, data_cpu = graph_cpu(data_frame, steps, hardware_info[host]['num_cpus'])

        graph.get_figure().savefig(
            os.path.join(outdir, '{}_cpu.png'.format(host)),
            bbox_inches='tight', pad_inches=0.25)
        pylab.close()

        ifaces = set([series.split('_')[0]
                      for series in data_frame.keys()
                      if series.startswith(('eth', 'ib'))])

        if verbose:
            print('Generating network graphs for {}...'.format(host))
        graph, data_net_bytes = graph_net_bytes(data_frame, steps, ifaces)
        graph.get_figure().savefig(
            os.path.join(outdir, '{}_net_bytes.png'.format(host)),
            bbox_inches='tight', pad_inches=0.25)
        pylab.close()

        graph, data_net_pkts = graph_net_pkts(data_frame, steps, ifaces)
        graph.get_figure().savefig(
            os.path.join(outdir, '{}_net_pkts.png'.format(host)),
            bbox_inches='tight', pad_inches=0.25)
        pylab.close()

        if verbose:
            print('Generating memory graph for {}...'.format(host))
        graph, data_mem = graph_memory(data_frame, steps, hardware_info[host]["memory"])
        graph.get_figure().savefig(
            os.path.join(outdir, '{}_memory.png'.format(host)),
            bbox_inches='tight', pad_inches=0.25)
        pylab.close()

        if verbose:
            print('Generating storage I/O graph for {}...'.format(host))
        drives = set([
            series.split('_')[0]
            for series in data_frame.keys()
            if series.startswith(('sd', 'vd', 'hd', 'xvd'))
        ])
        graph, data_disk = graph_disk_io(data_frame, steps, drives)
        graph.get_figure().savefig(
            os.path.join(outdir, '{}_disk_io.png'.format(host)),
            bbox_inches='tight', pad_inches=0.25)
        pylab.close()

        print('Serializing output to pickle object for node {}...'.format(host))
        # "Clean" dataframes ready to be plotted
        collectl_info[host] = { "hardware": hardware_info,
                                "steps": steps, "cpu": data_cpu, "mem": data_mem,
                                "disk": data_disk, "net_bytes": data_net_bytes,
                                "net_pkts": data_net_pkts
                              }
    return collectl_info

def serialize_plot_data(collectl_info, pre_graph_info, outdir, fname="collectl_info.pickle.gz"):
        # Useful to regenerate and slice graphs quickly and/or inspect locally
        collectl_pickle = os.path.join(outdir, fname)
        print("Saving plot pickle file with all hosts on: {}".format(collectl_pickle))
        with gzip.open(collectl_pickle, "wb") as f:
            pickle.dump((collectl_info, pre_graph_info), f)

def add_subparser(subparsers):
    parser = subparsers.add_parser(
        "graph",
        help=("Generate system graphs (CPU/memory/network/disk I/O "
              "consumption) from bcbio runs"))
    parser.add_argument(
        "log",
        help="Local path to bcbio log file written by the run.")
    parser.add_argument(
        "-o", "--outdir", default="monitoring/graphs",
        help="Directory to write graphs to.")
    parser.add_argument(
        "-r", "--rawdir", default="monitoring/collectl", required=True,
        help="Directory to put raw collectl data files.")
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Emit verbose output")

    return parser
