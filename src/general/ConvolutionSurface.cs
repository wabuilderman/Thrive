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
    public float GetIntegralAtPoint(Vector3 segA, Vector3 segB, Vector3 point, float sigma)
    {
        segA.X /= sigma;
        segA.Y /= sigma;
        segA.Z /= sigma;

        segB.X /= sigma;
        segB.Y /= sigma;
        segB.Z /= sigma;

        point.X /= sigma;
        point.Y /= sigma;
        point.Z /= sigma;

        int i = 3;
        int tau = 1;
        int tauDelta = 0;

        float integral = 0;
        for (int k = 0; k < i; k++)
        {
            float convolution = Convolution(point, segA, segB, i, k);
            integral += MathUtils.GetBinomialCoefficient(i, k)
                * Mathf.Pow(tauDelta, k) * Mathf.Pow(tau, i - k - 1) * convolution;
        }

        return integral / NormalizationFactor(i, sigma);
    }

    private float NormalizationFactor(int i, float sigma)
    {
        if (i == 2)
        {
            return sigma * sigma * Mathf.Pi;
        }
        else if (i == 3)
        {
            return sigma * sigma * sigma * 2;
        }
        else
        {
            return sigma * sigma * (i - 3) / (i - 2) * NormalizationFactor(i - 2, sigma);
        }
    }

    private float Convolution(Vector3 point, Vector3 segA, Vector3 segB, int i, int k)
    {
        float distAnB = Vector3.Distance(segA, segB);
        float distAnP = Vector3.Distance(segA, point);
        float distBnP = Vector3.Distance(segB, point);
        float dist2AnB = Vector3.DistanceSquared(segA, segB);
        float dist2AnP = Vector3.DistanceSquared(segA, point);
        float dist2BnP = Vector3.DistanceSquared(segB, point);
        Vector3 vectAnB = Vector3.Subtract(segA, segB);
        Vector3 vectAnP = Vector3.Subtract(segA, point);
        Vector3 vectBnP = Vector3.Subtract(segB, point);
        Vector3 vectBnA = Vector3.Subtract(segB, segA);

        float delta = dist2AnB * dist2AnP - Mathf.Pow(Vector3.Dot(vectAnB, vectAnP), 2);

        if (k == 0)
        {
            if (i == 1)
            {
                return Mathf.Log((distAnB * distBnP + Vector3.Dot(vectBnA, vectBnP))
                    / (distAnB * distAnP - Vector3.Dot(vectAnB, vectAnP)));
            }
            else if (i == 2)
            {
                return Mathf.Atan(Vector3.Dot(vectBnA, vectBnP / Mathf.Sqrt(delta))
                    + Mathf.Atan(Vector3.Dot(vectAnB, vectAnP / Mathf.Sqrt(delta)))
                    * distAnB / Mathf.Sqrt(delta));
            }
            else
            {
                return distAnB / (i - 2) / delta * ((i - 3) * dist2AnB * Convolution(point, segA, segB, i - 2, 0)
                    + Vector3.Dot(vectBnA, vectBnP) / Mathf.Pow(distBnP, i - 2)
                    + Vector3.Dot(vectAnB, vectAnP) / Mathf.Pow(distAnP, i - 2));
            }
        }
        else if (k == 1)
        {
            if (i == 2)
            {
                return Vector3.Dot(vectAnB, vectAnP) / dist2AnB * Convolution(point, segA, segB, 2, 0)
                    + Mathf.Log(dist2BnP / distAnP) / distAnB;
            }
            else
            {
                return Vector3.Dot(vectAnB, vectAnP) / dist2AnB * Convolution(point, segA, segB, i, i - 2)
                    + (Mathf.Pow(distBnP, 2 - i) - Mathf.Pow(distAnP, 2 - i)) / distAnB / (2 - i);
            }
        }
        else if (k == i - 1)
        {
            return Vector3.Dot(vectAnB, vectAnP) / dist2AnB * Convolution(point, segA, segB, i - 2, 0)
                + Convolution(point, segA, segB, i - 2, i - 3) / dist2AnB
                + 1 / ((2 - i) * distAnB * Mathf.Pow(distBnP, i - 2));
        }
        else
        {
            return ((i - 2 * k) * Vector3.Dot(vectAnB, vectAnP)) / ((i - k - 1) * dist2AnB)
                * Convolution(point, segA, segB, i, k - 1)
                + ((k - 1) * dist2AnP) / ((i - k - 1) * dist2AnB) * Convolution(point, segA, segB, i, k - 1)
                - 1 / ((distAnB * (i - k - 1)) * Mathf.Pow(distBnP, i - 2));
        }
    }
}
