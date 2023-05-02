using System;
using System.Collections.Generic;
using System.Linq;
using Godot;
using Supercluster.KDTree;

/// <summary>
///   Math related utility functions for Thrive
/// </summary>
public static class MathUtils
{
    public const float EPSILON = 0.00000001f;
    public const float DEGREES_TO_RADIANS = Mathf.Pi / 180;
    public const double FULL_CIRCLE = Math.PI * 2;
    public const float RIGHT_ANGLE = Mathf.Pi / 2;

    // standard euclidean distance
    private static readonly Func<float[], float[], double> l2Norm = (p, q) =>
    {
        double dist = 0;
        for (int i = 0; i < p.Length; i++)
        {
            dist += (p[i] - q[i]) * (p[i] - q[i]);
        }

        return dist;
    };

    public static T Clamp<T>(this T val, T min, T max)
        where T : IComparable<T>
    {
        if (val.CompareTo(min) < 0)
        {
            return min;
        }

        if (val.CompareTo(max) > 0)
        {
            return max;
        }

        return val;
    }

    public static double Sigmoid(double x)
    {
        return 1 / (1 + Math.Exp(-x));
    }

    /// <summary>
    ///   Creates a rotation for an organelle. This is used by the editor, but PlacedOrganelle uses RotateY as this
    ///   didn't work there for some reason.
    /// </summary>
    /// <remarks>
    ///   <para>
    ///     Rotation is now the number of 60 degree rotations
    ///   </para>
    /// </remarks>
    public static Quat CreateRotationForOrganelle(float rotation)
    {
        return new Quat(new Vector3(0, -1, 0), rotation * 60 * DEGREES_TO_RADIANS);
    }

    /// <summary>
    ///   This still takes the angle in degrees as this is used from
    ///   places that calculate the angle in degrees.
    /// </summary>
    public static Quat CreateRotationForExternal(float angle)
    {
        return new Quat(new Vector3(0, 1, 0), 180 * DEGREES_TO_RADIANS) *
            new Quat(new Vector3(0, 1, 0), angle * DEGREES_TO_RADIANS);
    }

    /// <summary>
    ///   Rotation for the pilus physics cone
    /// </summary>
    public static Quat CreateRotationForPhysicsOrganelle(float angle)
    {
        return new Quat(new Vector3(-1, 0, 0), 90 * DEGREES_TO_RADIANS) *
            new Quat(new Vector3(0, 0, -1), (180 - angle) * DEGREES_TO_RADIANS);
    }

    /// <summary>
    ///   Returns a Lerped value, and snaps to the target value if current and target
    ///   value is approximately equal by the specified tolerance value.
    /// </summary>
    public static float Lerp(float from, float to, float weight, float tolerance = Mathf.Epsilon)
    {
        if (Mathf.IsEqualApprox(from, to, tolerance))
            return to;

        return Mathf.Lerp(from, to, weight);
    }

    /// <summary>
    ///   Standard modulo for negative values in C# produces negative results.
    ///   This function returns modulo values between 0 and mod-1.
    /// </summary>
    /// <returns>The positive modulo</returns>
    public static int PositiveModulo(this int val, int mod)
    {
        int result = val % mod;
        return (result < 0) ? result + mod : result;
    }

    /// <summary>
    /// This function gets the total number of unique combinations based upon n and k.
    /// n is the total number of items.
    /// k is the size of the group.
    /// Total number of unique combinations = n! / ( k! (n - k)! ).
    /// This function is less efficient, but is more likely to not overflow when n and k are large.
    /// Taken from:  http://blog.plover.com/math/choose.html,
    /// https://stackoverflow.com/questions/12983731/algorithm-for-calculating-binomial-coefficient
    /// </summary>
    public static long GetBinomialCoefficient(long n, long k)
    {
        long r = 1;
        long d;
        if (k > n)
            return 0;

        for (d = 1; d <= k; d++)
        {
            r *= n--;
            r /= d;
        }

        return r;
    }

    public static (double Average, double StandardDeviation) CalculateAverageAndStandardDeviation(
        this IEnumerable<int> enumerable)
    {
        int count = 0;
        double sum = 0;
        double sumOfSquares = 0;

        foreach (var value in enumerable)
        {
            ++count;
            sum += value;
            sumOfSquares += value * value;
        }

        if (count == 0)
            throw new InvalidOperationException("Sequence contains no elements");

        double average = sum / count;
        double standardDeviation = Math.Sqrt(sumOfSquares / count - average * average);
        return (average, standardDeviation);
    }

    public static (double Average, double StandardDeviation) CalculateAverageAndStandardDeviation(
        this IEnumerable<double> enumerable)
    {
        int count = 0;
        double sum = 0;
        double sumOfSquares = 0;

        foreach (var value in enumerable)
        {
            ++count;
            sum += value;
            sumOfSquares += value * value;
        }

        if (count == 0)
            throw new InvalidOperationException("Sequence contains no elements");

        double average = sum / count;
        double standardDeviation = Math.Sqrt(sumOfSquares / count - average * average);
        return (average, standardDeviation);
    }

    /// <summary>
    ///   How far in a given direction a set of points can reach from a reference point
    /// </summary>
    /// <remarks>
    ///   <para>
    ///     This can be useful for finding the edge of a cell from a list of its organelle positions, for instance
    ///   </para>
    /// </remarks>
    public static float GetMaximumDistanceInDirection(Vector3 direction, Vector3 referencePoint,
        IEnumerable<Vector3> listOfPoints)
    {
        float distance = 0.0f;

        foreach (var point in listOfPoints)
        {
            if (point == referencePoint)
                continue;

            var difference = point - referencePoint;

            float angle = difference.AngleTo(direction);

            if (angle >= RIGHT_ANGLE)
                continue;

            // Get the length of the part of the vector that's parallel to the direction
            float directionalLength = difference.Length() * Mathf.Cos(angle);

            if (directionalLength > distance)
                distance = directionalLength;
        }

        return distance;
    }

    // find nearest neighbours using a kd tree
    private static IEnumerable<Vector3> GetNearestNeighbours(List<Vector3> sites)
    {
        int siteCount = sites.Count;

        var data = new List<float[]>();

        for (int i = 0; i < siteCount; i++)
        {
            data.Add(new float[] { sites[i].x, sites[i].y, sites[i].z });
        }

        float[][] treeData = data.ToArray();
        var treeNodes = sites.Select(p => p.ToString()).ToArray();
        var tree = new KDTree<float, string>(3, treeData, treeNodes, l2Norm);

        // gotta work on this one a lot
        throw new NotImplementedException();
    }

    // checks whether ray intersects triangle
    public static bool Intersect(Vector3 p1, Vector3 p2, Vector3 p3, Ray ray)
    {
        // Vectors from p1 to p2/p3 (edges)
        Vector3 edgeA, edgeB;

        Vector3 p, q, t;
        float determinant, invDeterminant, u, v;

        edgeA = p2 - p1;
        edgeB = p3 - p1;

        p = ray.Direction.Cross(edgeB);

        determinant = edgeA.Dot(p);

        // if determinant is near zero, ray lies in plane of triangle otherwise not
        if (determinant > -EPSILON && determinant < EPSILON)
            return false;

        invDeterminant = 1.0f / determinant;

        // calculate distance from p1 to ray origin
        t = ray.Origin - p1;

        // Calculate u parameter
        u = t.Dot(p) * invDeterminant;

        // Check for ray hit
        if (u is < 0 or > 1)
            return false;

        // Prepare to test v parameter
        q = t.Cross(edgeA);

        // Calculate v parameter
        v = ray.Direction.Dot(q) * invDeterminant;

        // Check for ray hit
        if (v < 0 || u + v > 1)
            return false;

        // ray does intersect
        if (edgeB.Dot(q) * invDeterminant > EPSILON)
            return true;

        // No hit at all
        return false;
    }

    public struct Ray
    {
        public Vector3 Origin;
        public Vector3 Direction;

        public Ray(Vector3 origin, Vector3 dir)
        {
            Origin = origin;
            Direction = dir;
        }
    }
}
