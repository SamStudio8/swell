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
    scheme_fh.seek(0)
    for line in scheme_fh:
        ref, start, end, tile, pool = line.strip().split()
        scheme, tile, side = tile.split("_", 2)
        if tiles[tile][0] != -1 and tiles[tile][1] != -1:
            l_tiles.append((scheme, tile, tiles[tile]))

    return l_tiles

def swell_from_depth(depth_path, tiles, genomes):
    depth_fh = open(depth_path)
    cursor = 0
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

    for t_i, (scheme_name, tile_num, tile) in enumerate(tiles):
        len_win = len(tile_dat[t_i])
        mean_cov = np.mean(tile_dat[t_i])
        median_cov = np.median(tile_dat[t_i])
        print(depth_path, tile_num, tile[0], tile[1], scheme_name, mean_cov, median_cov, len_win)

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
    group.add_argument("--bam")
    group.add_argument("--depth")
    parser.add_argument("--ref", required=True, nargs='+')
    parser.add_argument("--bed", required=False)

    args = parser.parse_args()

    if args.bed:
        tiles = load_scheme(args.bed)
    else:
        tiles = {}

    if args.bam:
        swell_from_bam(args.bam, tiles, args.ref[0])
    elif args.depth:
        swell_from_depth(args.depth, tiles, args.ref)

if __name__ == "__main__":
    main()
