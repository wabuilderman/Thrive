using System;
using System.Collections;
using System.Numerics;
using Godot;
using Array = System.Array;
using Vector3 = System.Numerics.Vector3;

public class ConvolutionSurface
{
    public static float GetIntegralAtPoint(Vector3 segA, Vector3 segB, Vector3 point, float sigma)
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

    private static float NormalizationFactor(int i, float sigma)
    {
        if (i == 2)
        {
            return sigma * sigma * Mathf.Pi;
        }

        if (i == 3)
        {
            return sigma * sigma * sigma * 2;
        }

        return sigma * sigma * (i - 3) / (i - 2) * NormalizationFactor(i - 2, sigma);
    }

    /// <summary>
    /// <para>E.Hubert, M-P.Cani 'Convolution Surfaces based on Polygonal Curve Skeletons'</para>
    /// HAL Open Science (https://hal.science) [Internet] 30 October 2009. id: inria-00429358
    /// <para>Available from: https://inria.hal.science/inria-00429358v1/document</para>
    /// </summary>
    /// <param name="point">a point in space to test against the surface; P</param>
    /// <param name="segA">Start of line segment associated with surface's skeleton; A</param>
    /// <param name="segB">End of line segment associated with surface's skeleton; B</param>
    /// <param name="i">Order of the kernel</param>
    /// <param name="k">For recurrence (See HAL Id: inria-00429358, section 3.2: Recurrences)</param>
    /// <remarks>
    ///    This function ultimately boils down to integrating:
    ///        \integral_0^1 (1 / (at^2-2bt+c)) ^ (1/2) dx<para/>
    ///
    ///    a= |AB|^2
    ///    <para>b= Vector(AB) * Vector(AP)</para>
    ///    c= |AP|^2
    ///    <para>t= time</para>
    ///
    ///    a-2b+c = |BP|^2
    ///    <para>a-b= Vector(BA) * Vector(BP)</para>
    ///    delta= |AB|^2 * |AP|^2 - (Vector(AB) * Vector(AP))^2
    /// </remarks>
    private static float Convolution(Vector3 point, Vector3 segA, Vector3 segB, int i, int k)
    {
        float distAnB = Vector3.Distance(segA, segB);
        float distAnP = Vector3.Distance(segA, point);
        float distBnP = Vector3.Distance(segB, point);
        float dist2AnB = Vector3.DistanceSquared(segA, segB);
        float dist2AnP = Vector3.DistanceSquared(segA, point);
        float dist2BnP = Vector3.DistanceSquared(segB, point);
        Vector3 vecAnB = Vector3.Subtract(segA, segB);
        Vector3 vecAnP = Vector3.Subtract(segA, point);
        Vector3 vecBnP = Vector3.Subtract(segB, point);
        Vector3 vecBnA = Vector3.Subtract(segB, segA);

        float delta = dist2AnB * dist2AnP - Mathf.Pow(Vector3.Dot(vecAnB, vecAnP), 2);

        if (k == 0)
        {
            if (i == 1)
            {
                return Mathf.Log((distAnB * distBnP + Vector3.Dot(vecBnA, vecBnP))
                    / (distAnB * distAnP - Vector3.Dot(vecAnB, vecAnP)));
            }

            if (i == 2)
            {
                return Mathf.Atan(Vector3.Dot(vecBnA, vecBnP / Mathf.Sqrt(delta))
                    + Mathf.Atan(Vector3.Dot(vecAnB, vecAnP / Mathf.Sqrt(delta)))
                    * distAnB / Mathf.Sqrt(delta));
            }

            return distAnB / (i - 2) / delta * ((i - 3) * dist2AnB * Convolution(point, segA, segB, i - 2, 0)
                + Vector3.Dot(vecBnA, vecBnP) / Mathf.Pow(distBnP, i - 2)
                + Vector3.Dot(vecAnB, vecAnP) / Mathf.Pow(distAnP, i - 2));
        }

        if (k == 1)
        {
            if (i == 2)
            {
                return Vector3.Dot(vecAnB, vecAnP) / dist2AnB * Convolution(point, segA, segB, 2, 0)
                    + Mathf.Log(dist2BnP / distAnP) / distAnB;
            }

            return Vector3.Dot(vecAnB, vecAnP) / dist2AnB * Convolution(point, segA, segB, i, i - 2)
                + (Mathf.Pow(distBnP, 2 - i) - Mathf.Pow(distAnP, 2 - i)) / distAnB / (2 - i);
        }

        if (k == i - 1)
        {
            return Vector3.Dot(vecAnB, vecAnP) / dist2AnB * Convolution(point, segA, segB, i - 2, 0)
                + Convolution(point, segA, segB, i - 2, i - 3) / dist2AnB
                + 1 / ((2 - i) * distAnB * Mathf.Pow(distBnP, i - 2));
        }

        return ((i - 2 * k) * Vector3.Dot(vecAnB, vecAnP)) / ((i - k - 1) * dist2AnB)
            * Convolution(point, segA, segB, i, k - 1)
            + ((k - 1) * dist2AnP) / ((i - k - 1) * dist2AnB) * Convolution(point, segA, segB, i, k - 1)
            - 1 / ((distAnB * (i - k - 1)) * Mathf.Pow(distBnP, i - 2));
    }
}
