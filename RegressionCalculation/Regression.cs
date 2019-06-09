using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq.Expressions;
using System.Text;
using static RegressionCalculation.Regression;
using static System.Math;

namespace RegressionCalculation
{
    /* TODO List
     * 1 add .FromParameters(A,B[,C]) (Values[])
     * 2 research Calculate and Find/Fill
     * 3 int NumOfParams
     * 4 ... etc
     */
    public abstract class Regression
    {
        #region Fields
        protected List<double> X, Y;
        protected int length;
        protected bool calculated;
        protected double r, a, b;
        #endregion
        #region Properties
        public bool Calculated => calculated;
        public double? R => calculated ? new double?(r) : null;
        public double? A => calculated ? new double?(a) : null;
        public double? B => calculated ? new double?(b) : null;
        public int Count => length;
        abstract protected string function { get; }
        public string Function => calculated ? function : "{not calculated}";
        #endregion
        #region ctors
        protected Regression()
        {
            X = new List<double>();
            Y = new List<double>();
        }
        protected Regression(double[] X, double[] Y)
        {
            length = X.Length;
            if (length != Y.Length) throw new Exception("Length");
            this.X = new List<double>(X);
            this.Y = new List<double>(Y);
        }
        #endregion
        #region Add-reset
        public void Add(double x, double y)
        {
            calculated = false;
            X.Add(x); Y.Add(y); length++;
        }
        public void AddRange(IEnumerable<double> x, IEnumerable<double> y)
        {
            calculated = false;
            X.AddRange(x); Y.AddRange(y);
            if ((length = X.Count) != Y.Count) throw new Exception("Length");
        }
        public void Reset() { X.Clear(); Y.Clear(); calculated = false; length = 0; }
        #endregion
        abstract public double GetY(double x);
        abstract public double GetX(double y);
        abstract protected void calculate();
        public Regression Calculate()
        {
            if (!calculated)
            {
                //ADDED FOR NOT THE FAKUP ON NAN/INF NUMBERS
                for (int i = length - 1; i >= 0; i--)
                    if (IsNaN(Y[i]))
                    {
                        Y.RemoveAt(i);
                        X.RemoveAt(i);
                        length--;
                    }
                //END
                calculate(); calculated = true;
                double res = 0, tot = 0, avg = 0;
                for (int i = 0; i < length; i++) avg += Y[i];
                avg /= length;
                for (int i = 0; i < length; i++)
                {
                    res += sqr(GetY(X[i]) - Y[i]);
                    tot += sqr(Y[i] - avg);
                }
                r = Sqrt(1 - res / tot);
            }
            a = Round(a, round);
            b = Round(b, round);
            r = Round(r, round);
            return this;
        }

        public virtual double[] GetValues() =>
            calculated ? new[] { r, a, b } : Calculate().GetValues();

        protected static double sqr(double v) => v * v;
        protected static bool IsNaN(double v) => double.IsInfinity(v) || double.IsNaN(v);
        protected static string[] Aformats = { "{0}", "", "+{0}" };
        protected static string[] Bformats = { "{1}{2}", "", "+{1}{2}" };
        protected static NumberFormatInfo nfi = new NumberFormatInfo() { NumberDecimalSeparator = "." };
        protected static int round = 10;

        static public int RoundDigit
        {
            get { return round; }
            set
            {
                if (value > 10 || value < -1) return;
                round = value;
            }
        }

        public override string ToString()
        {
            string str;
            if (!calculated) str = "{not calculated}";
            else
            {
                var tmp = GetValues();
                str = "{ R=" + tmp[0].ToString(nfi);
                for (int i = 1, n = tmp.Length; i < n; i++)
                    str += ", " + tmp[i].ToString(nfi);
                str += " }";
            }
            return GetType().Name + str;
        }
        internal void fromSource(List<double> x, List<double> y, int len)
        {
            length = len;
            X = x; Y = y;
            calculated = false;
        }
        public void FromSource(List<double> X, List<double> Y)
        {
            if (X.Count != Y.Count) throw new Exception("Length");
            fromSource(X, Y, X.Count);
        }

        public Regression FromParameters(double A, double B)
        {
            a = A; b = B; Reset(); calculated = true; r = 1;
            return this;
        }

        static B[] change<A, B>(A[] arr, Func<A, B> f)
        {
            int n = arr.GetLength(0);
            var ans = new B[n];
            for (int i = 0; i < n; i++)
                ans[i] = f(arr[i]);
            return ans;
        }
        static B[,] change<A, B>(A[,] arr, Func<A, B> f)
        {
            int n = arr.GetLength(0), m = arr.GetLength(1);
            var ans = new B[n, m];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                    ans[i, j] = f(arr[i, j]);
            return ans;
        }

        public virtual Func<double, double> GetFunction() => null;

        /// <summary>
        /// a + b * x
        /// </summary>
        public class Linear : Regression
        {
            public Linear() : base() { }
            public Linear(double[] X, double[] Y) : base(X, Y) { }

            public override double GetY(double x) => a + (b * x);
            public override double GetX(double y) => (y - a) / b;

            protected override void calculate()
            {
                double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, x, y;
                for (int i = 0; i < length; i++)
                {
                    x = X[i]; y = Y[i];
                    sumX += x; sumY += y;
                    sumXY += x * y;
                    sumX2 += x * x;
                }
                b = (length * sumXY - sumX * sumY) / (length * sumX2 - sumX * sumX);
                a = (sumY - b * sumX) / length;
            }
            protected override string function
                => string.Format(nfi, Aformats[Sign(a) + 1] + Bformats[Sign(b) + 1], a, b, "*x");

            public override Func<double, double> GetFunction()
            {
                var param = Expression.Parameter(typeof(double), "x");
                Expression<Func<double, double>> func = null;
                double A = a, B = b;
                func = x => A + B * x;
                func = Expression.Lambda<Func<double, double>>(Expression.Add(Expression.Constant(a), Expression.Multiply(Expression.Constant(b), param)), param);
                return func.Compile();
            }

            public double FindY(double x) => a + (b * x);
        }
        /// <summary>
        /// a + b * ln x
        /// </summary>
        public class Logarithmic : Regression
        {
            public Logarithmic(double[] X, double[] Y) : base(X, Y) { }
            public Logarithmic() : base() { }

            public override double GetY(double x) => a + b * Log(x);
            public override double GetX(double y) => Exp((y - a) / b);

            protected override void calculate()
            {
                double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, x, y;
                for (int i = 0; i < length; i++)
                {
                    x = Log(X[i]); y = Y[i];
                    sumX += x; sumY += y;
                    sumXY += x * y;
                    sumX2 += x * x;
                }
                b = (length * sumXY - sumX * sumY) / (length * sumX2 - sumX * sumX);
                a = (sumY - b * sumX) / length;
            }
            protected override string function
            => string.Format(nfi, Aformats[Sign(a) + 1] + Bformats[Sign(b) + 1], a, b, "*ln(x)");
        }
        /// <summary>
        /// a * e ^ (b * x)
        /// </summary>
        public class Exponential : Regression
        {
            public Exponential() : base() { }
            public Exponential(double[] X, double[] Y) : base(X, Y) { }

            public override double GetY(double x) => a * Exp(b * x);
            public override double GetX(double y) => Log(y / a) / b;

            protected override void calculate()
            {
                double x, y, sumX = 0, sumLnY = 0, sumXLnY = 0, sumX2 = 0;
                for (int i = 0; i < length; i++)
                {
                    x = X[i]; y = Y[i];
                    sumX += x;
                    sumLnY += Log(y);
                    sumX2 += x * x;
                    sumXLnY += x * Log(y);
                }
                b = (sumX * sumLnY - length * sumXLnY) / (sqr(sumX) - length * sumX2);
                a = Exp((sumLnY - b * sumX) / length);
            }
            protected override string function =>
                string.Format(nfi, "{0}*e^({1}*x)", a, b);
        }
        /// <summary>
        /// a * b ^ x
        /// </summary>
        public class TheExponential : Regression
        {
            public TheExponential() : base() { }
            public TheExponential(double[] X, double[] Y) : base(X, Y) { }

            public override double GetY(double x) => a * Pow(b, x);
            public override double GetX(double y) => Log(y / a) / Log(b);

            protected override void calculate()
            {
                double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, x, y;
                for (int i = 0; i < length; i++)
                {
                    x = X[i]; y = Log(Y[i]);
                    sumX += x; sumY += y;
                    sumXY += x * y;
                    sumX2 += x * x;
                }
                b = (length * sumXY - sumX * sumY) / (length * sumX2 - sumX * sumX);
                a = Exp((sumY - b * sumX) / length);
                b = Exp(b);
            }
            protected override string function
                => string.Format(nfi, "{0}*{1}^x", a, b);
        }
        /// <summary>
        /// a * x ^ b
        /// </summary>
        public class Power : Regression
        {
            public Power() : base() { }
            public Power(double[] X, double[] Y) : base(X, Y) { }

            public override double GetY(double x) => a * Pow(x, b);
            public override double GetX(double y) => Pow(y / a, 1 / b);

            protected override void calculate()
            {
                double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, x, y;
                for (int i = 0; i < length; i++)
                {
                    x = Log(X[i]); y = Log(Y[i]);
                    sumX += x; sumY += y;
                    sumXY += x * y;
                    sumX2 += x * x;
                }
                b = (length * sumXY - sumX * sumY) / (length * sumX2 - sumX * sumX);
                a = Exp((sumY - b * sumX) / length);
            }
            protected override string function
                => string.Format(nfi, "{0}*x^{1}", a, b);
        }
        /// <summary>
        /// a + b / x
        /// </summary>
        public class Inverse : Regression
        {
            public Inverse() : base() { }
            public Inverse(double[] X, double[] Y) : base(X, Y) { }

            public override double GetY(double x) => a + (b / x);
            public override double GetX(double y) => b / (y - a);

            protected override void calculate()
            {
                double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, x, y;
                for (int i = 0; i < length; i++)
                {
                    x = 1 / X[i]; y = Y[i];
                    sumX += x; sumY += y;
                    sumXY += x * y;
                    sumX2 += x * x;
                }
                b = (length * sumXY - sumX * sumY) / (length * sumX2 - sumX * sumX);
                a = (sumY - b * sumX) / length;
            }
            protected override string function
                => string.Format(nfi, Aformats[Sign(a) + 1] + Bformats[Sign(b) + 1], a, b, "/x");
        }
        /// <summary>
        /// a * x ^ 2 + b * x + c
        /// </summary>
        public class Quadratic : Regression
        {
            protected double c;
            public double? C => calculated ? new double?(c) : null;

            public override double GetY(double x) => a * x * x + b * x + c;
            public override double GetX(double y) => a == 0 ? (y - c) / b : (Sqrt(-4 * a * c + 4 * a * y + (b * b)) - b) / 2 / a;

            protected override void calculate()
            {
                double[] sumXpow = new double[5];
                sumXpow[0] = length;
                double x, y, t, sumY = 0, sumXY = 0, sumX2Y = 0;
                for (int i = 0; i < length; i++)
                {
                    x = X[i]; y = Y[i]; t = 1;
                    for (int pow = 1; pow <= 4; pow++)
                        sumXpow[pow] += (t *= x);
                    sumY += y;
                    sumXY += y * x;
                    sumX2Y += y * x * x;
                }
                double[,] sums = new double[3, 3];
                for (int ix = 0; ix < 3; ix++)
                    for (int iy = 0; iy < 3; iy++)
                        sums[ix, iy] = sumXpow[4 - ix - iy];
                Matrix values;
                {
                    var m = new Matrix(sums);
                    var a = m.Adjugate();
                    a = a.Transparent();
                    var d = m.Determinant();
                    m = a / d;
                    var n = Matrix.Vertical(new double[] { sumX2Y, sumXY, sumY });
                    values = m * n;
                }
                a = Round(values[0, 0], 10);
                b = Round(values[1, 0], 10);
                c = Round(values[2, 0], 10);
            }
            protected override string function
                => string.Format(nfi, "{0}*x^2+{1}*x+{2}", a, b, c);

            public override double[] GetValues() => new double[] { r, a, b, c };
        }
        /// <summary>
        /// a + b * sin(x)
        /// </summary>
        public class Sinus : Regression
        {
            public Sinus() : base() { }
            public Sinus(double[] X, double[] Y) : base(X, Y) { }

            public override double GetY(double x) => a + b * Sin(x);
            public override double GetX(double y) => Asin((y - a) / b);

            protected override void calculate()
            {
                double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, x, y;
                for (int i = 0; i < length; i++)
                {
                    x = Sin(X[i]); y = Y[i];
                    sumX += x; sumY += y;
                    sumXY += x * y;
                    sumX2 += x * x;
                }
                b = (length * sumXY - sumX * sumY) / (length * sumX2 - sumX * sumX);
                a = (sumY - b * sumX) / length;
            }
            protected override string function
                => string.Format(nfi, Aformats[Sign(a) + 1] + Bformats[Sign(b) + 1], a, b, "sin(x)");
        }
        /// <summary>
        /// a + b * cos(x)
        /// </summary>
        public class Cosinus : Regression
        {
            public Cosinus() : base() { }
            public Cosinus(double[] X, double[] Y) : base(X, Y) { }

            public override double GetY(double x) => a + b * Cos(x);
            public override double GetX(double y) => Acos((y - a) / b);

            protected override void calculate()
            {
                double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, x, y;
                for (int i = 0; i < length; i++)
                {
                    x = Cos(X[i]); y = Y[i];
                    sumX += x; sumY += y;
                    sumXY += x * y;
                    sumX2 += x * x;
                }
                b = (length * sumXY - sumX * sumY) / (length * sumX2 - sumX * sumX);
                a = (sumY - b * sumX) / length;
            }
            protected override string function
                => string.Format(nfi, Aformats[Sign(a) + 1] + Bformats[Sign(b) + 1], a, b, "cos(x)");
        }
        /// <summary>
        /// a + bx + cx2 + .. +zxn
        /// </summary>
        public class Polynomial : Regression
        {
            protected int m;
            protected double[] values;
            public int M => m;

            public Polynomial(int M) : base() { m = M; values = new double[m + 1]; }
            public Polynomial(double[] X, double[] Y, int M) : base(X, Y) { m = M; values = new double[m + 1]; }

            public override double GetY(double x)
            {
                double res = 0;
                for (int i = 0; i < m; i++)
                    res += values[i] * Pow(x, m - i);
                res += values[m];
                return res;
            }
            public override double GetX(double y) { throw new NotImplementedException(); }

            protected override void calculate()
            {
                int len = 2 * m - 2;
                double[]
                    sumXpow = new double[len + 1],
                    sumYpow = new double[m];
                sumXpow[0] = length;
                double x, y, t;
                for (int i = 0; i < length; i++)
                {
                    x = X[i]; y = Y[i]; t = 1;
                    for (int pow = 1; pow <= len; pow++)
                        sumXpow[pow] += (t *= x);
                    t = 1;
                    for (int pow = m - 1; pow >= 0; pow--)
                    {
                        sumYpow[pow] += y * t;
                        t *= x;
                    }
                }

                double[,] sums = new double[m, m];
                for (int ix = 0; ix < m; ix++)
                    for (int iy = 0; iy < m; iy++)
                        sums[ix, iy] = sumXpow[len - ix - iy];
#if DecMatrix
                var temp = new DecMatrix(change(sums, v => (decimal)v)).Invert() * DecMatrix.Vertical(change(sumYpow, v => (decimal)v));

                for (int i = 0; i < m; i++)
                    values[i + 1] = (double)Round(temp[i, 0], round);
#else
                var temp = new Matrix(sums).Invert() * Matrix.Vertical(sumYpow);

                for (int i = 0; i < m; i++)
                    values[i + 1] = Round(temp[i, 0], round);
#endif
                a = values[1];
                b = values[2];
            }
            protected override string function
            {
                get
                {
                    var sb = new StringBuilder();
                    for (int i = 1; i < m; i++)
                        sb.AppendFormat(Bformats[Sign(values[i]) + 1], null, values[i], "*x^" + (m - i));
                    return sb.AppendFormat(Aformats[Sign(values[m]) + 1], values[m]).ToString();
                }
            }
            public override double[] GetValues()
            {
                values[0] = r;
                return values;
            }
        }
        /// <summary>
        /// a + b/x + c/x2 + .. +z/xn
        /// </summary>
        public class InvertPolynomial : Regression
        {
            protected int m;
            protected double[] values;
            public int M => m;

            public InvertPolynomial(int M) : base() { m = M; values = new double[m + 1]; }
            public InvertPolynomial(double[] X, double[] Y, int M) : base(X, Y) { m = M; values = new double[m + 1]; }

            public override double GetY(double x)
            {
                double res = 0;
                for (int i = 0; i < m; i++)
                    res += values[i] / Pow(x, m - i);
                res += values[m];
                return res;
            }
            public override double GetX(double y)
            {
                throw new NotImplementedException();
            }

            protected override void calculate()
            {
                int len = 2 * m - 2;
                double[]
                    sumXpow = new double[len + 1],
                    sumYpow = new double[m];
                sumXpow[0] = length;
                double x, y, t;
                for (int i = 0; i < length; i++)
                {
                    x = 1 / X[i]; y = Y[i]; t = 1;
                    for (int pow = 1; pow <= len; pow++)
                        sumXpow[pow] += (t *= x);
                    t = 1;
                    for (int pow = m - 1; pow >= 0; pow--)
                    {
                        sumYpow[pow] += y * t;
                        t *= x;
                    }
                }
                double[,] sums = new double[m, m];
                for (int ix = 0; ix < m; ix++)
                    for (int iy = 0; iy < m; iy++)
                        sums[ix, iy] = sumXpow[len - ix - iy];
                var temp = new Matrix(sums).Invert() * Matrix.Vertical(sumYpow);
                for (int i = 0; i < m; i++)
                    values[i + 1] = Round(temp[i, 0], round);

                a = values[1];
                b = values[2];
            }
            protected override string function
            {
                get
                {
                    var sb = new StringBuilder();
                    for (int i = 1; i < m; i++)
                        sb.AppendFormat(Bformats[Sign(values[i]) + 1], null, values[i], "/x^" + (m - i));
                    return sb.AppendFormat(Aformats[Sign(values[m]) + 1], values[m]).ToString();
                }
            }
            public override double[] GetValues()
            {
                values[0] = r;
                return values;
            }
        }
    }
    public class FindRegression
    {
        const int length = 11;
        static Regression[] regressions = new Regression[length]
            {
                new Linear(),
                new Logarithmic(),
                new Exponential(),
                new TheExponential(),
                new Power(),
                new Inverse(),
                new Quadratic(),
                new Sinus(),
                new Cosinus(),
                new Polynomial(5),
                new InvertPolynomial(5)
            };
        static void fill()
        {
            //return;
            regressions = new Regression[length]
            {
                new Linear(),
                new Logarithmic(),
                new Exponential(),
                new TheExponential(),
                new Power(),
                new Inverse(),
                new Quadratic(),
                new Sinus(),
                new Cosinus(),
                new Polynomial(5),
                new InvertPolynomial(5)
            };
        }
        public static Regression[] Find(double[] X, double[] Y)
        {
            fill();
            int len, ind = -1, cnt = 0; double max = -1;
            if ((len = X.Length) != Y.Length) throw new Exception("Length");
            List<double> x = new List<double>(X), y = new List<double>(Y);
            for (int i = 0; i < length; i++)
            {
                var reg = regressions[i];
                reg.fromSource(x, y, len);
                reg.Calculate();
                if (reg.R > max) { max = reg.R.Value; ind = i; }
                if (reg.R == 1) cnt++;
            }
            if (cnt == 0) return ind >= 0 ? new Regression[] { regressions[ind] } : null;
            Regression[] tmp = new Regression[cnt];
            ind = 0;
            for (int i = 0; i < length; i++)
                if (regressions[i].R == 1)
                    tmp[ind++] = regressions[i];
            return tmp;
        }
        public static Regression[] Find(double[,] values)
        {
            int N = values.GetLength(0), M = values.GetLength(1);
            if (N != 2 && M != 2) throw new Exception("Length");
            double[] X, Y;
            if (N != 2)
            {
                X = new double[N];
                Y = new double[N];
                for (int i = 0; i < N; i++)
                {
                    X[i] = values[i, 0];
                    Y[i] = values[i, 1];
                }
            }
            else
            {
                X = new double[M];
                Y = new double[M];
                for (int i = 0; i < M; i++)
                {
                    X[i] = values[0, i];
                    Y[i] = values[1, i];
                }
            }
            return Find(X, Y);
        }
    }
}