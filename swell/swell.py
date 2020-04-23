import argparse
import sys
import pysam
import numpy as np

def load_scheme(scheme_bed):
    tiles = {}
    scheme_fh = open(scheme_bed)

    for line in scheme_fh:
        ref, start, end, tile, pool = line.strip().split()
        scheme, tile, side = tile.split("_", 2)

        if tile not in tiles:
            tiles[tile] = [-1, -1]

        if "LEFT" in side.upper():
            if tiles[tile][0] == -1:
                tiles[tile][0] = int(end)
            elif tiles[tile][0] < int(end):
                # If using multiple products for one tile, take the smallest window
                tiles[tile][0] = int(end)

        elif "RIGHT" in side.upper():
            if tiles[tile][1] == -1:
                tiles[tile][1] = int(start)
            elif tiles[tile][1] > int(start):
                tiles[tile][1] = int(start)

    l_tiles = []
    tiles_seen = set([])
    scheme_fh.seek(0)
    for line in scheme_fh:
        ref, start, end, tile, pool = line.strip().split()
        scheme, tile, side = tile.split("_", 2)
        tile_tup = (scheme, tile, tiles[tile])
        if tiles[tile][0] != -1 and tiles[tile][1] != -1 and tile_tup not in tiles_seen:
            l_tiles.append(tile_tup)
            tiles_seen.add(tile_tup)

    return l_tiles

def swell_from_fasta(fasta_path):
    num_seqs = 0
    num_bases = 0
    num_acgt = 0
    num_masked = 0
    num_invalid = 0

    prop_acgt = 0
    prop_masked = 0
    prop_invalid = 0

    if fasta_path:
        from . import readfq # thanks heng
        heng_iter = readfq.readfq(open(fasta_path))
        for name, seq, qual in heng_iter:
            num_seqs += 1
            for base in seq:
                num_bases += 1
                if base.upper() in ['A', 'C', 'G', 'T']:
                    num_acgt += 1
                elif base.upper() in ['N']:
                    num_masked += 1
                elif base.upper() in ['X', '-', '_', ' ']:
                    num_invalid += 1

        if num_bases > 0:
            prop_acgt = num_acgt / num_bases * 100.0
            prop_masked = num_masked / num_bases * 100.0
            prop_invalid = num_invalid / num_bases * 100.0
        else:
            prop_invalid = 100.0

    return ["fasta_path", "num_seqs", "num_bases", "pc_acgt", "pc_masked", "pc_invalid"], [fasta_path, num_seqs, num_bases, prop_acgt, prop_masked, prop_invalid]



def swell_from_depth(depth_path, tiles, genomes, thresholds):
    depth_fh = open(depth_path)

    threshold_counters = {
        threshold: 0 for threshold in thresholds
    }
    tile_threshold_counters = {
        threshold: 0 for threshold in thresholds
    }
    n_positions = 0
    avg_cov = 0

    cursor = 0
    if tiles:
        tile_starts = [t[2][0] for t in tiles] # dont use -1 for 1-pos depth files
        tile_ends = [t[2][1] for t in tiles]
        closest_cursor = min(tile_starts)

        stat_tiles = [0 for t in tiles]
        tile_dat = [[] for t in tiles]

    for line in depth_fh:
        ref, pos, cov = line.strip().split('\t')
        if sum([g in ref for g in genomes]) != 1:
            continue
        pos = int(pos)
        cov = int(cov)

        # Count positions above threshold
        for threshold in threshold_counters:
            if cov >= threshold:
                threshold_counters[threshold] += 1
        n_positions += 1
        avg_cov = avg_cov + (cov - avg_cov)/n_positions

        if tiles:
            # Check for new open tiles
            if pos >= closest_cursor:
                for t_i, t_start in enumerate(tile_starts):
                    if pos >= t_start and stat_tiles[t_i] == 0:
                        stat_tiles[t_i] = 1

                next_possible_min = []
                for t_i, t_start in enumerate(tile_starts):
                    if stat_tiles[t_i] == 0:
                        next_possible_min.append(t_start)
                try:
                    closest_cursor = min(next_possible_min)
                except ValueError:
                    closest_cursor = sys.maxsize

            # Handle open tiles
            for t_i, t_state in enumerate(stat_tiles):
                if t_state == 1:
                    tile_dat[t_i].append(cov)

                if tile_ends[t_i] <= pos:
                    stat_tiles[t_i] = -1

    tile_vector = []
    for t_i, (scheme_name, tile_num, tile) in enumerate(tiles):
        len_win = len(tile_dat[t_i])
        mean_cov = np.mean(tile_dat[t_i])
        median_cov = np.median(tile_dat[t_i])

        tile_vector.append(mean_cov)

        # Count tile means above threshold
        for threshold in threshold_counters:
            if mean_cov >= threshold:
                tile_threshold_counters[threshold] += 1

        #print(depth_path, tile_num, tile[0], tile[1], scheme_name, mean_cov, median_cov, len_win)


    if n_positions > 0:
        threshold_counts_prop = [threshold_counters[x]/n_positions * 100.0 for x in sorted(thresholds)]
    else:
        threshold_counts_prop = [0 for x in sorted(thresholds)]

    if len(tiles) > 0:
        tile_threshold_counts_prop = [tile_threshold_counters[x]/len(tiles) * 100.0 for x in sorted(thresholds)]
    else:
        tile_threshold_counts_prop = [0 for x in sorted(thresholds)]

    if len(tile_vector) == 0:
        tile_vector.append(0)

    return ["bam_path", "num_pos", "mean_cov"] + ["pc_pos_cov_gte%d" % x for x in sorted(thresholds)] + ["pc_tiles_meancov_gte%d" % x for x in sorted(thresholds)] + ["tile_vector"], [depth_path.replace(".depth", ""), n_positions, avg_cov] + threshold_counts_prop + tile_threshold_counts_prop + [",".join(["%.2f" % x for x in tile_vector])]

#def swell_from_bam(bam_path, tiles, genome):
#    bam = pysam.AlignmentFile(bam_path)
#
#    for (scheme_name, tile_num, tile) in tiles:
#        tile_cover = bam.count_coverage(genome, tile[0]-1, tile[1],
#                quality_threshold=0, read_callback="all")
#        flat_tile_cover = np.array(tile_cover).sum(axis=0)
#
#        mean_cov = np.mean(flat_tile_cover)
#        median_cov = np.median(flat_tile_cover)
#        print(bam_path, tile_num, tile[0], tile[1], scheme_name, mean_cov, median_cov)

def main():
    import argparse
    parser = argparse.ArgumentParser()
#    group = parser.add_mutually_exclusive_group(required=True)
#    group.add_argument("--bam")
    parser.add_argument("--depth", required=True)
    parser.add_argument("--ref", required=True, nargs='+')
    parser.add_argument("--thresholds", action='append', type=int, nargs='+', default=[1, 5, 10, 20, 50, 100, 200])
    parser.add_argument("--bed", required=False)
    parser.add_argument("--fasta", required=False)
    parser.add_argument("--dp", default=2, type=int, required=False)

    args = parser.parse_args()

    if args.bed:
        tiles = load_scheme(args.bed)
    else:
        tiles = {}

    #if args.bam:
    #    swell_from_bam(args.bam, tiles, args.ref[0])
    #elif args.depth:
    #    swell_from_depth(args.depth, tiles, args.ref)
    fields = []
    header = []

    header_, fields_ = swell_from_fasta(args.fasta)
    header.extend(header_)
    fields.extend(fields_)

    header_, fields_ = swell_from_depth(args.depth, tiles, args.ref, args.thresholds)
    header.extend(header_)
    fields.extend(fields_)

    print("\t".join(header))
    fields_s = [("%."+str(args.dp)+"f") % x if "float" in type(x).__name__ else str(x) for x in fields] # do not fucking @ me
    print("\t".join([str(x) for x in fields_s]))

if __name__ == "__main__":
    main()
