#!/usr/bin/env python
# Python implementation of Sebastien's example.

if __name__ == '__main__':
    import sys
    from string import atof

    try:
        with open(sys.argv[1]) as f:
            line = f.readline().split()
            t0, v0, t1, v1 = [atof(line[i]) for i in range(4)]
    except:
        sys.exit(1)

    t2 = int(3 - t0 - t1)
    print '%d %d %d %d' % (t2, v0, t1, v1)  # Neighbor 1.
    print '%d %d %d %d' % (t0, v0, t2, v1)  # Neighbor 2.

