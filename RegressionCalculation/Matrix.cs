using System;

public class Matrix
{
    [Serializable]
    public class MatrixException : Exception
    {
        public MatrixException() { }
    }

    protected double[,] matrix;
    protected int n, m;
    public int N => n;
    public int M => m;

    public Matrix(double[,] arr)
    {
        matrix = arr;
        n = arr.GetLength(0);
        m = arr.GetLength(1);
    }
    public Matrix(int N, int M)
    {
        matrix = new double[N, M];
        n = N; m = M;
    }
    protected Matrix(int N, int M, double[,] matrix)
    {
        n = N; m = M;
        this.matrix = matrix;
    }

    public double this[int x, int y]
    {
        get { return matrix[x, y]; }
        set { matrix[x, y] = value; }
    }

    public bool Square => n == m;
    public static bool Like(Matrix a, Matrix b) => a.n == b.n && a.m == b.m;

    public static Matrix operator +(Matrix a, Matrix b)
    {
        if (!Like(a, b)) throw new MatrixException();
        int n = a.n, m = a.m;
        var tmp = new double[n, m];
        for (int x = 0; x < n; x++)
            for (int y = 0; y < m; y++)
                tmp[x, y] = a.matrix[x, y] + b.matrix[x, y];
        return new Matrix(n, m, tmp);
    }
    public static Matrix operator -(Matrix a, Matrix b)
    {
        if (!Like(a, b)) throw new MatrixException();
        int n = a.n, m = a.m;
        var tmp = new double[n, m];
        for (int x = 0; x < n; x++)
            for (int y = 0; y < m; y++)
                tmp[x, y] = a.matrix[x, y] - b.matrix[x, y];
        return new Matrix(n, m, tmp);
    }
    public static Matrix operator *(Matrix a, Matrix b)
    {
        if (a.m != b.n) throw new MatrixException();
        double[,] tmp = new double[a.n, b.m];
        double sum;
        for (int x = a.n - 1; x >= 0; --x)
        {
            for (int y = b.m - 1; y >= 0; --y)
            {
                sum = 0;
                for (int i = a.m - 1; i >= 0; --i)
                {
                    sum += a.matrix[x, i] * b.matrix[i, y];
                }
                tmp[x, y] = sum;
            }
        }
        return new Matrix(tmp);
    }
    public static Matrix operator *(Matrix a, double v)
    {
        int n = a.n, m = a.m;
        var tmp = new double[n, m];
        for (int x = 0; x < n; x++)
            for (int y = 0; y < m; y++)
                tmp[x, y] = a.matrix[x, y] * v;
        return new Matrix(n, m, tmp);
    }
    public static Matrix operator /(Matrix a, double v)
    {
        int n = a.n, m = a.m;
        var tmp = new double[n, m];
        for (int x = 0; x < n; x++)
            for (int y = 0; y < m; y++)
                tmp[x, y] = a.matrix[x, y] / v;
        return new Matrix(n, m, tmp);
    }

    public Matrix AddMinor(int cutX, int cutY)
    {
        var tmp = new double[n - 1, m - 1];
        for (int x = 0, xx = 0; xx < n - 1; x++, xx++)
        {
            if (x == cutX) x++;
            for (int y = 0, yy = 0; yy < m - 1; y++, yy++)
            {
                if (y == cutY) y++;
                tmp[xx, yy] = matrix[x, y];
            }
        }
        return new Matrix(n - 1, m - 1, tmp);
    }
    public double Minor(int x, int y) => AddMinor(x, y).Determinant();
    public double Addition(int x, int y) => Minor(x, y) * (((x + y + 1) & 1) * 2 - 1);
    public double Determinant()
    {
        if (n != m || n <= 0) throw new MatrixException();
        if (n == 1) return matrix[0, 0];
        if (n == 2) return matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0];
        if (n == 3)
            return
            matrix[0, 0] * matrix[1, 1] * matrix[2, 2] - matrix[0, 0] * matrix[1, 2] * matrix[2, 1] +
            matrix[0, 1] * matrix[1, 2] * matrix[2, 0] - matrix[0, 1] * matrix[1, 0] * matrix[2, 2] +
            matrix[0, 2] * matrix[1, 0] * matrix[2, 1] - matrix[0, 2] * matrix[1, 1] * matrix[2, 0];
        double ans = 0;
        for (int x = 0, i = 1; x < n; x++, i = -i)
            ans += i * matrix[x, 0] * Minor(x, 0);
        return ans;
    }
    public Matrix Transparent()
    {
        double[,] tmp = new double[m, n];
        for (int x = 0; x < n; x++)
            for (int y = 0; y < m; y++)
                tmp[y, x] = matrix[x, y];
        return new Matrix(m, n, tmp);

    }
    public Matrix Adjugate()
    {
        double[,] tmp = new double[n, m];
        for (int x = 0; x < n; x++)
            for (int y = 0; y < m; y++)
                tmp[x, y] = Addition(x, y);
        return new Matrix(n, m, tmp);
    }
    public Matrix Invert()
    {
        return Adjugate().Transparent() / Determinant();
    }
    public double permanent()
    {
        if (n == 1) return matrix[0, 0];
        if (n == 2) return matrix[0, 0] * matrix[1, 1] + matrix[0, 1] * matrix[1, 0];

        double ans = 0;
        for (int x = 0; x < n; x++)
            ans += matrix[x, 0] * AddMinor(x, 0).permanent();
        return ans;
    }
    public double Permanent()
    {
        if (!Square) throw new MatrixException();
        return permanent();
    }

    public static Matrix Horizontal(double[] v)
    {
        int m = v.Length;
        var tmp = new double[1, m];
        for (int i = 0; i < m; i++)
            tmp[0, i] = v[i];
        return new Matrix(1, m, tmp);
    }
    public static Matrix Vertical(double[] v)
    {
        int n = v.Length;
        var tmp = new double[n, 1];
        for (int i = 0; i < n; i++)
            tmp[i, 0] = v[i];
        return new Matrix(n, 1, tmp);
    }
}
public class DecMatrix
{
    protected decimal[,] matrix;
    protected int n, m;
    public int N => n;
    public int M => m;

    public DecMatrix(decimal[,] arr)
    {
        matrix = arr;
        n = arr.GetLength(0);
        m = arr.GetLength(1);
    }
    public DecMatrix(int N, int M)
    {
        matrix = new decimal[N, M];
        n = N; m = M;
    }
    protected DecMatrix(int N, int M, decimal[,] matrix)
    {
        n = N; m = M;
        this.matrix = matrix;
    }

    public decimal this[int x, int y]
    {
        get { return matrix[x, y]; }
        set { matrix[x, y] = value; }
    }

    public bool Square => n == m;
    public static bool Like(DecMatrix a, DecMatrix b) => a.n == b.n && a.m == b.m;

    public static DecMatrix operator +(DecMatrix a, DecMatrix b)
    {
        if (!Like(a, b)) throw new Exception();
        int n = a.n, m = a.m;
        var tmp = new decimal[n, m];
        for (int x = 0; x < n; x++)
            for (int y = 0; y < m; y++)
                tmp[x, y] = a.matrix[x, y] + b.matrix[x, y];
        return new DecMatrix(n, m, tmp);
    }
    public static DecMatrix operator -(DecMatrix a, DecMatrix b)
    {
        if (!Like(a, b)) throw new Exception();
        int n = a.n, m = a.m;
        var tmp = new decimal[n, m];
        for (int x = 0; x < n; x++)
            for (int y = 0; y < m; y++)
                tmp[x, y] = a.matrix[x, y] - b.matrix[x, y];
        return new DecMatrix(n, m, tmp);
    }
    public static DecMatrix operator *(DecMatrix a, DecMatrix b)
    {
        if (a.m != b.n) throw new Exception();
        decimal[,] tmp = new decimal[a.n, b.m];
        decimal sum;
        for (int x = 0; x < a.n; x++)
        {
            for (int y = 0; y < b.m; y++)
            {
                sum = 0;
                for (int i = 0; i < a.m; i++)
                {
                    sum += a.matrix[x, i] * b.matrix[i, y];
                }
                tmp[x, y] = sum;
            }
        }
        return new DecMatrix(tmp);
    }
    public static DecMatrix operator *(DecMatrix a, decimal v)
    {
        int n = a.n, m = a.m;
        var tmp = new decimal[n, m];
        for (int x = 0; x < n; x++)
            for (int y = 0; y < m; y++)
                tmp[x, y] = a.matrix[x, y] * v;
        return new DecMatrix(n, m, tmp);
    }
    public static DecMatrix operator /(DecMatrix a, decimal v)
    {
        int n = a.n, m = a.m;
        var tmp = new decimal[n, m];
        for (int x = 0; x < n; x++)
            for (int y = 0; y < m; y++)
                tmp[x, y] = a.matrix[x, y] / v;
        return new DecMatrix(n, m, tmp);
    }

    public DecMatrix AddMinor(int cutX, int cutY)
    {
        var tmp = new decimal[n - 1, m - 1];
        for (int x = 0, xx = 0; xx < n - 1; x++, xx++)
        {
            if (x == cutX) x++;
            for (int y = 0, yy = 0; yy < m - 1; y++, yy++)
            {
                if (y == cutY) y++;
                tmp[xx, yy] = matrix[x, y];
            }
        }
        return new DecMatrix(n - 1, m - 1, tmp);
    }
    public decimal Minor(int x, int y) => AddMinor(x, y).Determinant();
    public decimal Addition(int x, int y) => Minor(x, y) * (((x + y + 1) & 1) * 2 - 1);
    public decimal Determinant()
    {
        if (n != m || n <= 0) throw new Exception();
        if (n == 1) return matrix[0, 0];
        if (n == 2) return matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0];
        if (n == 3)
            return
            matrix[0, 0] * matrix[1, 1] * matrix[2, 2] - matrix[0, 0] * matrix[1, 2] * matrix[2, 1] +
            matrix[0, 1] * matrix[1, 2] * matrix[2, 0] - matrix[0, 1] * matrix[1, 0] * matrix[2, 2] +
            matrix[0, 2] * matrix[1, 0] * matrix[2, 1] - matrix[0, 2] * matrix[1, 1] * matrix[2, 0];
        decimal ans = 0;
        for (int x = 0, i = 1; x < n; x++, i = -i)
            ans += i * matrix[x, 0] * Minor(x, 0);
        return ans;
    }
    public DecMatrix Transparent()
    {
        decimal[,] tmp = new decimal[m, n];
        for (int x = 0; x < n; x++)
            for (int y = 0; y < m; y++)
                tmp[y, x] = matrix[x, y];
        return new DecMatrix(m, n, tmp);

    }
    public DecMatrix Adjugate()
    {
        decimal[,] tmp = new decimal[n, m];
        for (int x = 0; x < n; x++)
            for (int y = 0; y < m; y++)
                tmp[x, y] = Addition(x, y);
        return new DecMatrix(n, m, tmp);
    }
    public DecMatrix Invert()
    {
        return Adjugate().Transparent() / Determinant();
    }
    public decimal permanent()
    {
        if (n == 1) return matrix[0, 0];
        if (n == 2) return matrix[0, 0] * matrix[1, 1] + matrix[0, 1] * matrix[1, 0];

        decimal ans = 0;
        for (int x = 0; x < n; x++)
            ans += matrix[x, 0] * AddMinor(x, 0).permanent();
        return ans;
    }
    public decimal Permanent()
    {
        if (!Square) throw new Exception();
        return permanent();
    }

    public static DecMatrix Horizontal(decimal[] v)
    {
        int m = v.Length;
        var tmp = new decimal[1, m];
        for (int i = 0; i < m; i++)
            tmp[0, i] = v[i];
        return new DecMatrix(1, m, tmp);
    }
    public static DecMatrix Vertical(decimal[] v)
    {
        int n = v.Length;
        var tmp = new decimal[n, 1];
        for (int i = 0; i < n; i++)
            tmp[i, 0] = v[i];
        return new DecMatrix(n, 1, tmp);
    }
}
