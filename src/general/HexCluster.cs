using System;
using System.Collections.Generic;
using Godot;

public class HexCluster
{
    public List<Hex> MaxGridSpace(int maxHexes)
    {
        List<Hex> grid = new()
        {
            // this should be center of mass, ideally
            [0] = new Hex(0, 0),
        };

        int sliceArea = 0;
        int offset = 0;

        // offset is how big the side of a slice will be
        // we can use that to find the nth triangle number
        while (sliceArea < maxHexes / 6)
        {
            offset++;
            sliceArea = (offset * offset - offset) / 2;
        }

        int offsetQ = offset;

        // makes triangle of hexes like a pie slice of the grid
        for (int r = 1; r <= offset; r++)
        {
            for (int q = 1; q <= offsetQ; q++)
            {
                // 6 way symmetry
                grid.Add(new Hex(q, r));
                grid.Add(new Hex(-1 * r, r + q));
                grid.Add(new Hex(-1 * (r + q), q));
                grid.Add(new Hex(-1 * q, -1 * r));
                grid.Add(new Hex(r, -1 * (r + q)));
                grid.Add(new Hex(r + q, -1 * q));
            }
        }

        return grid;
    }
}
