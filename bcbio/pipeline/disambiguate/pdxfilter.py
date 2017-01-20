
import time
import re
import sys
import pysam
import logging
import os
import errno
from operator import itemgetter
from argparse import ArgumentParser, RawTextHelpFormatter
from collections import defaultdict

logging_formater = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger           = logging.getLogger()

logger.setLevel(logging.ERROR)

stream_handler   = logging.StreamHandler()

stream_handler.setFormatter(logging_formater)
logger.addHandler(stream_handler)

PATTERN_INTEGER = re.compile('([0-9]+)')


def mkdir(path):
    dirpath = os.path.dirname(path)
    try:
        os.makedirs(dirpath)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(dirpath):
            pass
        else:
            raise


def convert_to_digit(string):
    if string.isdigit():
        return int(string)
    return string


def convert_alphanumeric(string):
    return [convert_to_digit(substring)
            for substring in PATTERN_INTEGER.split(string)]


def greater_natural(name_1, name_2):
    """ SAMtools sorts queries by mixed alphanumeric sort rather than by lexicographic
        sort. We emulate this behaviour here.
    """

    return convert_alphanumeric(name_1) > convert_alphanumeric(name_2)


def greater_lexicographic(name_1, name_2):
    """ Sambamba sorts query by lexiocographic sort.
    """

    return name_1 > name_2


class PDXFilter(object):

    def __init__(self, human_sam_path, explant_sam_path,
                 output_human_sam_path, output_explant_sam_path,
                 output_ambiguous_sam_path, summary_path=None,
                 hard_filter=True, debug=False):

        self.input_human_sam_path   = human_sam_path
        self.input_explant_sam_path = explant_sam_path

        self.human_sam_file   = pysam.AlignmentFile(self.input_human_sam_path, 'r')
        self.explant_sam_file = pysam.AlignmentFile(self.input_explant_sam_path, 'r')

        # If alignments should be filtered hard or soft
        # A hard filter decides the fates of all alignments of a read pair based
        # on the highest-scoring alignment.
        # Instead, a soft filter splits all human alignments of a read based on the
        # best-scoring alignment of that read againt the explant reference.
        # All human alignments scoring higher than the best explant alignment
        # go to human while all other alignments go to ambiguous.
        self.hard_filter = hard_filter

        # Ordering scheme of the SAM/BAM files
        self.greater = None

        # Set output paths
        self.output_human_sam_path      = output_human_sam_path
        self.output_explant_sam_path    = output_explant_sam_path
        self.output_ambiguous_sam_path  = output_ambiguous_sam_path
        self.summary_path               = summary_path

        # Make output directories
        mkdir(self.output_human_sam_path)
        mkdir(self.output_explant_sam_path)
        mkdir(self.output_ambiguous_sam_path)

        if summary_path:
            mkdir(summary_path)

        output_human_mode       = output_human_sam_path.endswith('.bam') and 'wb' or 'wh'
        output_explant_mode     = output_explant_sam_path.endswith('.bam') and 'wb' or 'wh'
        output_ambiguous_mode   = output_ambiguous_sam_path.endswith('.bam') and 'wb' or 'wh'

        self.output_human_sam_file     = pysam.AlignmentFile(self.output_human_sam_path,
                                                             mode=output_human_mode,
                                                             template=self.human_sam_file)

        self.output_explant_sam_file   = pysam.AlignmentFile(self.output_explant_sam_path,
                                                             mode=output_explant_mode,
                                                             template=self.human_sam_file)

        self.output_ambiguous_sam_file = pysam.AlignmentFile(self.output_ambiguous_sam_path,
                                                             mode=output_ambiguous_mode,
                                                             template=self.human_sam_file)

        # Summary information
        self.human_bundles_count     = 0
        self.explant_bundles_count   = 0
        self.ambiguous_bundles_count = 0

        self.debug = debug

        if self.debug:
            logger.setLevel(logging.DEBUG)
        else:
            sys.tracebacklimit = 0

    def set_sort_order(self, alignments):

        natural_sort       = True
        lexicographic_sort = True
        current_name       = None

        for alignment in alignments:

            # Set first query name
            if not current_name:

                current_name = alignment.qname
                continue

            elif current_name != alignment.qname:

                # Check for natural sort order
                if not greater_natural(alignment.qname, current_name):
                    natural_sort = False

                # Check for lexicographic sort order
                if not greater_lexicographic(alignment.qname, current_name):
                    lexicographic_sort = False

                if not natural_sort and not lexicographic_sort:
                    break

                current_name = alignment.qname

        # We prefer natural sort since use of samtools is more common
        if natural_sort:
            logging.info('Disambiguation: Assuming natural (mixed) sort order of both input SAM/BAM files.')
            self.greater = greater_natural

        elif lexicographic_sort:
            logging.info('Disambiguation: Assuming lexicographic sort order of both input SAM/BAM files.')
            self.greater = greater_lexicographic

        # If the ordering does not correspond to any known order
        else:
            raise RuntimeError('Disambiguation: SAM/BAM file does not appear to be ordered by query name. Using human alignment {} and explant alignment {}.'.format(self.input_human_sam_path, self.input_explant_sam_path))

    def iter_sam_file(self, sam_file, check_first=10000):

        alignments = []
        num_alignments = 0
        try:
		    for i, alignment in enumerate(sam_file.fetch(until_eof=True,
		                                                 multiple_iterators=False)):

		        num_alignments += 1

		        # If we know the sort order, we continue
		        if self.greater:
		            yield(alignment)
		            continue

		        # Else w cache the first N alignments
		        elif i < check_first:
		            alignments.append(alignment)
		            continue

		        # Or infer the sort order baed on the cache
		        else:
		            self.set_sort_order(alignments)

		            # We can yield the stored alignments now that the sort order is inferred
		            for cached_alignment in alignments:
		                yield cached_alignment
		            alignments = []
		            yield alignment

		    # If we had less than check_first alignments in the file, we have to determine
		    # sort order based on that
		    if not self.greater:
		        self.set_sort_order(alignments)
		        for cached_alignment in alignments:
		            yield cached_alignment

		    logging.info('Disambiguation: a SAM file resulted in {} alignments'.format(num_alignments))
		    
        except:
            logging.error('Disambiguation: Error with human alignment {} and/or explant alignment {}.'.format(self.input_human_sam_path, self.input_explant_sam_path))
            raise

    def get_query_bundles(self, sam_file):
        """ Generatator that bundles multiple (multimapping) alignment of the
            same query read. Does not seperate paired ends within the bundles,
            i.e. bundles contain paired ends in an arbitrary order.
        """

        num_bundles    = 0
        num_alignments = 0
        bundle       = []
        current_name = None

        for alignment in self.iter_sam_file(sam_file):

            num_alignments += 1

            # Initialize
            if current_name is None:

                current_name = alignment.qname
                bundle       = [alignment]
                continue

            # Additional alignment for the current bundle
            elif alignment.qname == current_name:
                bundle.append(alignment)

            # A new bundle (alignments with another read name) begins
            else:

                # Fall back to lexicographic ordering
                if not self.greater(alignment.qname, current_name):

                    message = 'Disambiguation: SAM/BAM file does not appear to be ordered by'\
                              ' query name ({} is not greater {})'\
                              ' Using human alignment {} and explant alignment {}.'\
                              .format(alignment.qname,
                                      current_name,
                                      self.input_human_sam_path,
                                      self.input_explant_sam_path)

                    if self.greater is greater_natural and \
                            greater_lexicographic(alignment.qname, current_name):
                        self.greater = greater_lexicographic
                        logging.warn(message + ' Falling back to lexicographic ordering.')
                    else:
                        raise ValueError(message)

                if not bundle:
                    logging.warn('Disambiguation: Produced empty bundle based on current name {}, alignment_name {}, and alignment {}. Produced {} alignments and {} valid bundles before.'.format(current_name, alignment.qname, alignment, num_alignments, num_bundles))
                else:
                    num_bundles += 1

                yield (current_name, bundle)

                # Rewrite information with new bundle
                current_name = alignment.qname
                bundle       = [alignment]

        if not bundle:
            logging.warn('Disambiguation: Produced empty last bundle based on current name {}. Produced {} alignments and {} valid bundles before.'.format(current_name, num_alignments, num_bundles))
        else:
            num_bundles += 1

        if num_alignments == 0:
            logging.warn('Disambiguation: Did not receive any alignments from either human SAM {} or explant SAM {}'.format(self.input_human_sam_path, self.input_explant_sam_path))

        if num_bundles == 0:
            logging.warn('Disambiguation: Did not receive any alignment bundles from either human SAM {} or explant SAM {}'.format(self.input_human_sam_path, self.input_explant_sam_path))

        yield (current_name, bundle)

    def get_paired_query_bundles(self, sam_file):
        """ Generatator that bundles multiple (multimapping) alignment of the
            same query read. Seperates paired ends within the bundles,
            i.e. bundles contain tuples of paired ends. Note that paired ends
            are not necessarily located in succeeding lines even if the SAM
            file is sorted by query name; therefore, we have to match paired
            ends manually across the whole query bundle.
        """

        for name, bundle in self.get_query_bundles(sam_file):

            # Do sanity check on query names
            query_names = set([alignment.qname for alignment in bundle])
            valid       = len(query_names) == 1 and name == list(query_names)[0]

            if not valid:
                message = 'Alignments incorrectly bundled. No or several query names ({}) in bundle {} for human SAM {}, explant SAM {}'.format(query_names, bundle, self.input_human_sam_path, self.input_explant_sam_path)
                raise IOError(message)

            alignment_keys = defaultdict(set)
            paired_bundle  = []

            for alignment in bundle:

                # For an unpaired alignment, we can directly store it as a bundle
                if not alignment.is_paired:

                    if alignment.is_read1:
                        paired_bundle.append((alignment, None))
                    else:
                        paired_bundle.append((None, alignment))

                # Paired alignments have to be reconstructed so that we have
                # both mates
                else:

                    self_mapped = not alignment.is_unmapped
                    mate_mapped = not alignment.mate_is_unmapped

                    # Summarize the current alignment
                    # We pad integer values by one to work with 0 values
                    # and retain the shortened logical expression
                    self_key = (alignment.is_read1,
                                self_mapped and
                                (alignment.reference_id + 1) or None,
                                self_mapped and
                                (alignment.reference_start + 1) or None,
                                mate_mapped and
                                (alignment.next_reference_id + 1) or None,
                                mate_mapped and
                                (alignment.next_reference_start + 1) or None)

                    # Describe how the mate of the current alignment
                    # would look like.
                    mate_key = (not alignment.is_read1,
                                mate_mapped and
                                (alignment.next_reference_id + 1) or None,
                                mate_mapped and
                                (alignment.next_reference_start + 1) or None,
                                self_mapped and
                                (alignment.reference_id + 1) or None,
                                self_mapped and
                                (alignment.reference_start + 1) or None)

                    # Try to recover stored mates for the current alignment
                    if mate_key not in alignment_keys:

                        # If we cannot retrieve the mate for a segment we have to
                        # store the alignment (not its mate) for later
                        alignment_keys[self_key].add(alignment)

                    # Else we have mates and are ready to pair them up with the alignment
                    else:

                        for mate_alignment in alignment_keys[mate_key]:

                            # If succeeded, build a read pair
                            if mate_alignment.is_read1:
                                paired_bundle.append((mate_alignment, alignment))

                            else:
                                paired_bundle.append((alignment, mate_alignment))

            if not paired_bundle:
                message = 'Disambiguation: Encountered empty bundle of alignments in SAM/BAM'\
                          ' file at read name {}. Mate pairing dailed and this'\
                          ' read will be skipped from the analysis.'.format(name)
                logging.warn(message)
                continue

            yield (name, paired_bundle)

    def get_matched_query_bundles(self):

        human_queries   = iter(self.get_paired_query_bundles(self.human_sam_file))
        explant_queries = iter(self.get_paired_query_bundles(self.explant_sam_file))

        last_explant_name          = None
        last_explant_paired_bundle = None

        for human_name, human_paired_bundle in human_queries:

            # We have found the matching records from the over-stepped
            # record in the previous search
            if last_explant_name == human_name:
                yield (human_paired_bundle, last_explant_paired_bundle)
                continue

            # Seek to matching explant query record
            while True:

                # We arrived at the last explant record without finding a match
                try:
                    explant_name, explant_paired_bundle = next(explant_queries)

                except StopIteration:

                    yield (human_paired_bundle, None)
                    break

                # We found the matching records
                if explant_name == human_name:

                    yield (human_paired_bundle, explant_paired_bundle)
                    break

                # We sought one step too far without finding a matching record,
                # producing one over-stepped record
                elif self.greater(explant_name, human_name):

                    last_explant_name          = explant_name
                    last_explant_paired_bundle = explant_paired_bundle

                    # Also, we can directly return the human record since it
                    # does not have an explant counterpart
                    yield (human_paired_bundle, None)
                    break

                # We have not sought far enough across the explant reads
                else:
                    continue

    def score_bundle(self, bundle):

        scored_bundle = []
        for end_1, end_2 in bundle:

            # Single.end or paired-end unmapped
            if not end_1 or (end_1 and end_1.is_unmapped):
                AS_end_1 = 0
                NM_end_1 = 0.0
                nM_end_1 = 0.0

            # Paired-end mapped
            else:
                tags_1      = dict(end_1.tags)
                AS_end_1    = tags_1.get('AS', 0)
                NM_end_1    = 1.0 / (tags_1.get('NM', 0) + 1)
                nM_end_1    = 1.0 / (tags_1.get('nM', 0) + 1)

            # Single.end or paired-end unmapped
            if not end_2 or (end_2 and end_2.is_unmapped):
                AS_end_2 = 0
                NM_end_2 = 0.0
                nM_end_2 = 0.0

            # Paired-end mapped
            else:
                tags_2      = dict(end_2.tags)
                AS_end_2    = tags_2.get('AS', 0)
                NM_end_2    = 1.0 / (tags_2.get('NM', 0) + 1)
                nM_end_2    = 1.0 / (tags_2.get('nM', 0) + 1)

            # Accumulate scores into a tuple
            score = (AS_end_1 + AS_end_2,
                     NM_end_1 + NM_end_2 +
                     nM_end_1 + nM_end_2)

            record = (score, (end_1, end_2))
            scored_bundle.append(record)

        sorted_bundle = sorted(scored_bundle, reverse=True, key=itemgetter(0))

        return sorted_bundle

    def write_bundle(self, bundle, sam_file):

        written = set()

        for end_1, end_2 in bundle:

            if end_1 and not end_1 in written:
                sam_file.write(end_1)
                written.add(end_1)

            if end_2 and not end_2 in written:
                sam_file.write(end_2)
                written.add(end_2)

    def set_corrected_tags(self, bundle):

        # If no alignment survives filtering, we do not have to regenerate tags
        number_alignments = len(bundle)

        # Now we can regenerate the tags to specify the number of surviving
        # alignments
        for i, (end_1, end_2) in enumerate(bundle):

            if end_1:

                if not end_1.is_unmapped:
                    end_1.setTag('NH', number_alignments)
                    end_1.setTag('HI', i+1)
                else:
                    end_1.setTag('NH', 0)
                    end_1.setTag('HI', 1)

            if end_2:

                if not end_2.is_unmapped:
                    end_2.setTag('NH', number_alignments)
                    end_2.setTag('HI', i+1)
                else:
                    end_2.setTag('NH', 0)
                    end_2.setTag('HI', 1)

    def disambiguate(self):

        for human_paired_bundle, explant_paired_bundle\
                in self.get_matched_query_bundles():

            # Resulting sorted read bundles that are written to SAM files
            human_filtered_bundle     = []
            explant_filtered_bundle   = []
            ambiguous_filtered_bundle = []

            # If the human bundle does not have an explant correspondence,
            # we can instantly write the segments to file without retagging
            if explant_paired_bundle is None:
                human_filtered_bundle = human_paired_bundle

            # Else we have to score and disambiguate
            else:

                # Score paired alignments within the bundles
                human_scored_bundle     = self.score_bundle(human_paired_bundle)
                explant_scored_bundle   = self.score_bundle(explant_paired_bundle)

                # Extract best scores for human and explant
                best_human_score        = human_scored_bundle[0][0]
                best_explant_score      = explant_scored_bundle[0][0]

                # If the best explant read pair is unmapped, do not filter based on it
                # but directly assign all human alignment to human
                if sum(best_explant_score) == 0:
                    human_filtered_bundle  = human_paired_bundle

                # Decide based on a hard filter, i.e. assign all human
                # alignments based on the the scores of the best human and
                # explant alignments
                elif self.hard_filter:

                    if best_human_score > best_explant_score:
                        human_filtered_bundle = human_paired_bundle

                    elif best_human_score == best_explant_score:
                        ambiguous_filtered_bundle = human_paired_bundle

                    else:
                        explant_filtered_bundle = human_paired_bundle

                # Decide based on a soft filter, i.e. split all human
                # alignments based on the the scores of the best explant alignment
                else:

                    for (human_score, (end_1, end_2)) in human_scored_bundle:

                        if human_score > best_explant_score:
                            human_filtered_bundle.append((end_1, end_2))

                        elif human_score == best_explant_score:
                            ambiguous_filtered_bundle.append((end_1, end_2))

                        else:
                            explant_filtered_bundle.append((end_1, end_2))

                    # Since we split the alignments we now have to retag
                    # them to specifiy the number of alignments for each
                    # read pair
                    if ambiguous_filtered_bundle or explant_filtered_bundle:
                        if human_filtered_bundle:
                            self.set_corrected_tags(human_filtered_bundle)
                        if ambiguous_filtered_bundle:
                            self.set_corrected_tags(ambiguous_filtered_bundle)
                        if explant_filtered_bundle:
                            self.set_corrected_tags(explant_filtered_bundle)

            yield human_filtered_bundle, ambiguous_filtered_bundle, explant_filtered_bundle

    def run(self):

        # Prepare debug output
        last_time       = None
        start_time      = None
        alignments_all  = 0
        alignments_last = 0

        for human_filtered_bundle,\
            ambiguous_filtered_bundle,\
                explant_filtered_bundle in self.disambiguate():

            # Now we can update the statistics based on how the paired bundle
            # was distributed.
            if human_filtered_bundle:
                self.human_bundles_count += len(human_filtered_bundle)

            if explant_filtered_bundle:
                self.explant_bundles_count += len(explant_filtered_bundle)

            if ambiguous_filtered_bundle:
                self.ambiguous_bundles_count += len(ambiguous_filtered_bundle)

            # Write alignments to the corresponding output files
            if human_filtered_bundle:
                self.write_bundle(human_filtered_bundle,
                                  self.output_human_sam_file)

            if ambiguous_filtered_bundle:
                self.write_bundle(ambiguous_filtered_bundle,
                                  self.output_ambiguous_sam_file)

            if explant_filtered_bundle:
                self.write_bundle(explant_filtered_bundle,
                                  self.output_explant_sam_file)

            # Handle debug counting
            if self.debug:

                number_alignments = (len(human_filtered_bundle) +
                                     len(ambiguous_filtered_bundle) +
                                     len(explant_filtered_bundle))

                alignments_all  += number_alignments
                alignments_last += number_alignments

                now             = time.time()

                if not start_time:
                    start_time = now
                    last_time  = now

                time_diff = (now - last_time)

                if time_diff > 60:

                    overal_time_diff = (now - start_time) / 3600.0
                    last_time_diff   = (now - last_time)  / 3600.0

                    logging.debug('Disambiguation: Runtime %5.1fh\tAnalyzed %5iM alignment pairs;\tAnalyzing at average rates of\t%5iM alignment pairs/h (based on overall throughput)\t%5iM alignment pairs/h (based on last minute)' % (
                                  overal_time_diff, alignments_all / 10**6, (alignments_all / overal_time_diff) / 10 ** 6, (alignments_last / last_time_diff) / 10 ** 6))

                    last_time       = now
                    alignments_last = 0

        self.output_human_sam_file.close()
        self.output_ambiguous_sam_file.close()
        self.output_explant_sam_file.close()

        if self.summary_path:

            with open(self.summary_path, 'w') as summary_file:

                header = ['sample', 'unique species A pairs',
                          'unique species B pairs', 'ambiguous pairs']
                summary_file.write('\t'.join(header) + '\n')

                sample_name = os.path.splitext(os.path.basename(self.input_human_sam_path))[0]

                content = [sample_name,
                           self.human_bundles_count,
                           self.explant_bundles_count,
                           self.ambiguous_bundles_count]

                summary_file.write('\t'.join(map(str, content)) + '\n')

if __name__ == '__main__':

    description = """ Filters explant reads from a humam SAM/BAM alignment using
                      a second alignment of the same reads against the explant species.

                      Input files:
                      human_sam_path: Path to input SAM/BAM file of reads mapped against human
                      explant_sam_path: Path to input SAM/BAM file of reads mapped against explant species

                      Output files:
                      output_human_sam_path: Output to SAM/BAM file (decided by suffix) of cleaned human reads
                      output_explant_sam_path: Output to SAM/BAM file (decided by suffix) of cleaned human reads
                      output_ambiguous_sam_path: Output to SAM/BAM file (decided by suffix) of ambiguous reads that were removed from the 'human' outputs
                      summary_path (optional): Output to text-delimited text file with filtering statistics
                    """

    parser = ArgumentParser(description=description,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument('--human_sam_path', help='Path to input SAM/BAM file of reads mapped against human')
    parser.add_argument('--explant_sam_path', help='Path to input SAM/BAM file of reads mapped against explant species.')
    parser.add_argument('--output_human_sam_path',
                        help='Output to SAM/BAM file (decided by suffix) of cleaned human reads')
    parser.add_argument('--output_explant_sam_path',
                        help='Output to SAM/BAM file (decided by suffix) of explant reads that were removed from the human outputs')
    parser.add_argument('--output_ambiguous_sam_path',
                        help='Output to SAM/BAM file (decided by suffix) of ambiguous reads that were removed from the human outputs')
    parser.add_argument('--soft_filter', action='store_true', default=False,
                        help='Enable soft filtering instead of hard filtering')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='Enable debug output')
    parser.add_argument('--summary_path', default='',
                        help='Output to text-delimited text file with filtering statistics')
    args = parser.parse_args()

    hard_filter = not args.soft_filter

    pdx_filter = PDXFilter(args.human_sam_path,
                           args.explant_sam_path,
                           args.output_human_sam_path,
                           args.output_explant_sam_path,
                           args.output_ambiguous_sam_path,
                           args.summary_path,
                           hard_filter=hard_filter,
                           debug=args.debug)
    pdx_filter.run()
