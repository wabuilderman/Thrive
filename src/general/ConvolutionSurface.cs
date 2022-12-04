using System;
using System.Collections;
using System.Numerics;
using Godot;
using Array = System.Array;
using Vector3 = System.Numerics.Vector3;

/// <summary>
///    point: a point in space to test against the surface; P
///    segA: Start of line segment associated with surface's skeleton; A
///    segB: End of line segment associated with surface's skeleton; B
///    i: Order of the kernel
///    k: For determining recurrences (See HAL Id: inria-00429358, section 3.2: Recurrences)
///
///    This function ultimately boils down to integrating:
///        \integral_0^1 (1 / (at^2-2bt+c)) ^ (1/2) dx
///
///    a= |AB|^2
///    b= Vector(AB) * Vector(AP)
///    c= |AP|^2
///    t= time
///
///    a-2b+c = |BP|^2
///    a-b= Vector(BA) * Vector(BP)
///    delta= |AB|^2 * |AP|^2 - (Vector(AB) * Vector(AP))^2
/// </summary>
public class ConvolutionSurface
{
    private float Convolution(Vector3 point, Vector3 segA, Vector3 segB, int i, int k)
    {
        float distAB = Vector3.Distance(segA, segB);
        float distAP = Vector3.Distance(segA, point);
        float distBP = Vector3.Distance(segB, point);
        float dist2AB = Vector3.DistanceSquared(segA, segB);
        float dist2AP = Vector3.DistanceSquared(segA, point);
        float dist2BP = Vector3.DistanceSquared(segB, point);
        Vector3 vectAB = Vector3.Subtract(segA, segB);
        Vector3 vectAP = Vector3.Subtract(segA, point);
        Vector3 vectBP = Vector3.Subtract(segB, point);
        Vector3 vectBA = Vector3.Subtract(segB, segA);

        float delta = dist2AB * dist2AP - Mathf.Pow(Vector3.Dot(vectAB, vectAP), 2);

        if (k == 0)
        {
            if (i == 1)
            {
                return Mathf.Log((distAB * distBP + Vector3.Dot(vectBA, vectBP))
                    / (distAB * distAP - Vector3.Dot(vectAB, vectAP)));
            }
            else if (i == 2)
            {
                return Mathf.Atan(Vector3.Dot(vectBA, vectBP / Mathf.Sqrt(delta))
                    + Mathf.Atan(Vector3.Dot(vectAB, vectAP / Mathf.Sqrt(delta)))
                    * distAB / Mathf.Sqrt(delta));
            }
            else
            {
                return distAB / (i - 2) / delta * ((i - 3) * dist2AB * Convolution(point, segA, segB, i - 2, 0)
                    + Vector3.Dot(vectBA, vectBP) / Mathf.Pow(distBP, i - 2)
                    + Vector3.Dot(vectAB, vectAP) / Mathf.Pow(distAP, i - 2));
            }
        }
        else if (k == 1)
        {
            if (i == 2)
            {
                return Vector3.Dot(vectAB, vectAP) / dist2AB * Convolution(point, segA, segB, 2, 0)
                    + Mathf.Log(dist2BP / distAP) / distAB;
            }
            else
            {
                return Vector3.Dot(vectAB, vectAP) / dist2AB * Convolution(point, segA, segB, i, i - 2)
                    + (Mathf.Pow(distBP, 2 - i) - Mathf.Pow(distAP, 2 - i)) / distAB / (2 - i);
            }
        }
        else if (k == i - 1)
        {
            return Vector3.Dot(vectAB, vectAP) / dist2AB * Convolution(point, segA, segB, i - 2, 0)
                + Convolution(point, segA, segB, i - 2, i - 3) / dist2AB
                + 1 / ((2 - i) * distAB * Mathf.Pow(distBP, i - 2));
        }
        else
        {
            return ((i - 2 * k) * Vector3.Dot(vectAB, vectAP)) / ((i - k - 1) * dist2AB)
                * Convolution(point, segA, segB, i, k - 1)
                + ((k - 1) * dist2AP) / ((i - k - 1) * dist2AB) * Convolution(point, segA, segB, i, k - 1)
                - 1 / ((distAB * (i - k - 1)) * Mathf.Pow(distBP, i - 2));
        }
    }
}
