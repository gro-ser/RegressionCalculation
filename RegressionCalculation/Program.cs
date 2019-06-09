using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Linq.Expressions;
using System.Drawing;

using static System.Math;
using static RegressionCalculation.Regression;
using Microsoft.Win32;

namespace RegressionCalculation
{
    class Program
    {
        private static int len = 50;
        static Random R = new Random();
        static double f(double x) => 5 * Pow(x, -2);
        static void ExprTest()
        {
            const int a = 5, b = 3;
            Expression<Func<int>> g;
            g = () => a + b;
            //g = Expression.Lambda<Func<int>>(Expression.Add(Expression.Constant(a), Expression.Constant(b)));
            var f = g.Compile();
            Console.WriteLine(f());
            //a = -2;
            Console.WriteLine(f());
        }
        static void forTest(float s)
        {
            float A, B, C;
            A = -2; B = 2; C = 1;
            for (var i = A * s; i <= B * s; i += s)
                Console.Write(i / s + " ");
            Console.WriteLine();
        }
        static void Main(string[] args)
        {
            int length = 5;
            void results(Func<double, double> f)
            {
                for (int x = 0; x < length; x++)
                    Console.Write("{0:+0;-0;0}\t", f(x-2));
            }
            var fun = new Linear();
            fun.FromParameters(3, 5);
            Func<double, double>[] funcs =
            { fun.FindY, fun.GetFunction1(), fun.GetFunction2(), fun.GetFunction3(), fun.GetFunction4() };
            Console.WriteLine("Function string:" + fun.Function);
            for (int i = 0; i < length; i++)
            {
                Console.Write($"func{i}:\t");
                results(funcs[i]);
                Console.WriteLine(funcs[i].Target);
            }
            fun.FromParameters(5, 3);
            Console.WriteLine();
            Console.WriteLine("Function string:" + fun.Function);
            for (int i = 0; i < length; i++)
            {
                Console.Write($"func{i}:\t");
                results(funcs[i]);
                Console.WriteLine(funcs[i].Method);
            }
            Console.ReadLine();
        }

        long map(long x, long in_min, long in_max, long out_min, long out_max)
        {
            return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
        }

        static public double[] interpolate(double[] arr, int len)
        {
            int cnt = arr.Length;
            mapper map = new mapper(0, cnt - 1, 1, len);
            double[] X = new double[cnt], Y = new double[cnt], ans = new double[len];
            for (int i = 0; i < cnt; i++)
            {
                X[i] = map.map(i);
                Y[i] = arr[i];
            }
            var t = FindRegression.Find(X, Y);
            Console.Error.WriteLine("LENGTH:" + t.Length);
            var f = t[0];
            Console.Error.WriteLine("FUNCTION:" + f.Function);
            for (int i = 0; i < len; i++)
            {
                ans[i] = f.GetY(i + 1);
            }
            return ans;
        }

        static Regression[] StrToFun(string str)
        => FindRegression.Find(str.Select((a, i) => (double)i).ToArray(),
            str.Select(x => (double)x).ToArray());

        static double f1(double x)
            => -0.189231112 * Pow(x, 4) + 4.3524441245 * Pow(x, 3) - 31.0566553555 * x * x + 69.8611217191 * x + 67.2776346379;
        static double f2(double x) => 0.5 * x * x + 0.3 * x + 49.3;

        static string FunToString(Func<double, double> f, int len = 15)
        {
            var sb = new StringBuilder();
            for (int i = 0; i < len; i++)
                sb.Append((char)f(i));
            return sb.ToString();
        }

        static void FFF(string str, int len = 100)
        {
            var f = StrToFun(str)[0];
            Console.WriteLine(f);
            Console.WriteLine(f.Function);
            Console.WriteLine(FunToString(f.GetY, len));
        }

        static Regression fontest(int low = 5, int hi = 50)
        {
            var font = new FontFamily("Consolas");
            double[] arr0 = new double[hi - low], arr1 = new double[hi - low]; ;
            for (int i = low; i < hi; i++)
            {
                arr0[i - low] = i;
                arr1[i - low] = new Font(font, i).Height;
            }
            return new Regression.Linear(arr0, arr1).Calculate();
        }
    }
    class mapper
    {
        double in_min, in_max, out_min, out_max;
        public mapper(double in_min, double in_max, double out_min, double out_max)
        {
            this.in_min = in_min;
            this.in_max = in_max;
            this.out_min = out_min;
            this.out_max = out_max;
        }
        public double map(double x) => (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
    }
}
