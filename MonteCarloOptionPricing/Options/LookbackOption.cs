using System;
using System.Diagnostics;

// This file is in collaboration with my classmate Zongxing Li. 

namespace MonteCarloOptionPricing.Options
{
    public class LookbackOption
    {
        static double phai(double x)
        {
            double c1 = 0.254829592;
            double c2 = -0.284496736;
            double c3 = 1.421413741;
            double c4 = -1.453152027;
            double c5 = 1.061405429;
            double p = 0.3275911;
            double t, y;

            int c_sign = 1;
            if (x < 0)
                c_sign = -1;
            x = Math.Abs(x) / Math.Sqrt(2.0);

            t = 1.0 / (1.0 + p * x);
            y = 1.0 - (((((c5 * t + c4) * t) + c3) * t + c2) * t + c1) * t * Math.Exp(-x * x);

            return 0.5 * (1.0 + c_sign * y);
        }

        public static (double option_value, double standard_error) option_origin(double S, double K, double T, double r, double sigma, long N, int steps, string call_put, string multi_thread)
        {
            double price = 0;
            double standard_error = 0;
            double dt = T / steps;

            double[] option_value = new double[N];     
            double[,] stock_price = new double[N, steps];    
            double[,] seq = new double[N, steps];
            double[] stock_price_max = new double[N];   
            double[] stock_price_min = new double[N];   

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < steps; j++)
                {
                    seq[i, j] = Random_Number.Box_Muller_Transform();
                }
            }

            if (multi_thread == "no")
            {
                for (int i = 0; i < N; i++)
                {
                    stock_price[i, 0] = S;
                    stock_price_max[i] = S;
                    stock_price_min[i] = S;

                    for (int j = 1; j < steps; j++)
                    {
                        stock_price[i, j] = stock_price[i, j - 1] * Math.Exp((r - 0.5 * sigma * sigma) * dt + sigma * Math.Sqrt(dt) * seq[i, j]);
                        stock_price_max[i] = Math.Max(stock_price_max[i], stock_price[i, j]);
                        stock_price_min[i] = Math.Min(stock_price_min[i], stock_price[i, j]);
                    }
                    if (call_put == "call")    
                    {
                        option_value[i] = Math.Exp(-r * T) * Math.Max(stock_price_max[i] - K, 0);
                    }
                    else                  
                    {
                        option_value[i] = Math.Exp(-r * T) * Math.Max(K - stock_price_min[i], 0);
                    }
                }
            }
            else
            {
                int coreCount = Environment.ProcessorCount;
                ParallelOptions options = new ParallelOptions();
                options.MaxDegreeOfParallelism = coreCount;
                Parallel.For(0, N, options, i =>
                {
                    stock_price[i, 0] = S;
                    stock_price_max[i] = S;
                    stock_price_min[i] = S;

                    for (int j = 1; j < steps; j++)
                    {
                        stock_price[i, j] = stock_price[i, j - 1] * Math.Exp((r - 0.5 * sigma * sigma) * dt + sigma * Math.Sqrt(dt) * seq[i, j]);
                        stock_price_max[i] = Math.Max(stock_price_max[i], stock_price[i, j]);
                        stock_price_min[i] = Math.Min(stock_price_min[i], stock_price[i, j]);
                    }
                    if (call_put == "call")    
                    {
                        option_value[i] = Math.Exp(-r * T) * Math.Max(stock_price_max[i] - K, 0);
                    }
                    else                  
                    {
                        option_value[i] = Math.Exp(-r * T) * Math.Max(K - stock_price_min[i], 0);
                    }
                });
            }

            price = option_value.Average();
            standard_error = Math.Sqrt(option_value.Select(x => Math.Pow(x - option_value.Average(), 2)).Sum() / (N - 1)) / Math.Sqrt(N);

            return (price, standard_error);    
        }


        public static (double delta, double gamma, double vega, double theta, double rho) option_origin_greeks(double S, double K, double T, double r, double sigma, long N, int steps, string call_put, string multi_thread)
        {
            double delta, delta1, delta2, gamma, vega, theta, rho;
            delta = (option_origin(S * 1.05, K, T, r, sigma, N, steps, call_put, multi_thread).option_value - option_origin(S * 0.95, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (2 * S * 0.05);

            delta1 = (option_origin(S * 1.05, K, T, r, sigma, N, steps, call_put, multi_thread).option_value - option_origin(S, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (S * 0.05);
            delta2 = (option_origin(S, K, T, r, sigma, N, steps, call_put, multi_thread).option_value - option_origin(S * 0.95, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (S * 0.05);
            gamma = 2 * (delta1 - delta2) / (2 * S * 0.05);

            vega = (option_origin(S, K, T, r, sigma * 1.05, N, steps, call_put, multi_thread).option_value - option_origin(S, K, T, r, sigma * 0.95, N, steps, call_put, multi_thread).option_value) / (sigma * 1.05 - sigma * 0.95);
            theta = (option_origin(S, K, T * 0.95, r, sigma, N, steps, call_put, multi_thread).option_value - option_origin(S, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (T * 0.05);
            rho = (option_origin(S, K, T, r * 1.05, sigma, N, steps, call_put, multi_thread).option_value - option_origin(S, K, T, r * 0.95, sigma, N, steps, call_put, multi_thread).option_value) / (r * 1.05 - r * 0.95);

            return (delta, gamma, vega, theta, rho);
        }

        static double option_mc_delta(double S, double K, double T, double r, double sigma, string call_put)
        {
            double detla;
            double d1 = (Math.Log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * Math.Sqrt(T));
            if (call_put == "call")
            {
                detla = phai(d1);
            }
            else
            {
                detla = phai(d1) - 1;
            }
            return detla;
        }

        public static (double option_value, double standard_error) option_cv(double S, double K, double T, double r, double sigma, long N, int steps, string call_put, string multi_thread)
        {
            double price = 0;
            double standard_error = 0;
            double dt = T / steps;

            double nudt = (r - 0.5 * sigma * sigma) * dt;
            double sigmadt = sigma * Math.Sqrt(dt);
            double erddt = Math.Exp(r * dt);
            double beta1 = -1;

            double[] option_value = new double[N];     
            double[,] stock_price = new double[N, steps];    
            double[,] seq = new double[N, steps];
            double[,] cv = new double[N, steps];
            double[] CT = new double[N];
            double[,] local_t = new double[N, steps];
            double[,] delta = new double[N, steps];
            double[] stock_price_max = new double[N];
            double[] stock_price_min = new double[N];

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < steps; j++)
                {
                    seq[i, j] = Random_Number.Box_Muller_Transform();
                }
            }

            if (multi_thread == "no")
            {
                for (int i = 0; i < N; i++)
                {
                    stock_price[i, 0] = S;
                    stock_price_max[i] = S;
                    stock_price_min[i] = S;

                    for (int j = 1; j < steps; j++)
                    {
                        local_t[i, j] = j * dt;
                        delta[i, j] = option_mc_delta(stock_price[i, j - 1], K, T - local_t[i, j - 1], r, sigma, call_put);
                        stock_price[i, j] = stock_price[i, j - 1] * Math.Exp(nudt + sigmadt * seq[i, j]);
                        cv[i, j] = delta[i, j] * (stock_price[i, j] - stock_price[i, j - 1] * erddt);

                        stock_price_max[i] = Math.Max(stock_price_max[i], stock_price[i, j]);
                        stock_price_min[i] = Math.Min(stock_price_min[i], stock_price[i, j]);
                    }

                    if (call_put == "call")    
                    {
                        CT[i] = Math.Max(0, stock_price_max[i] - K) + beta1 * cv[i, steps - 1];
                    }

                    else                  
                    {
                        CT[i] = Math.Max(0, K - stock_price_min[i]) + beta1 * cv[i, steps - 1];
                    }
                }
            }
            else
            {
                int coreCount = Environment.ProcessorCount;
                ParallelOptions options = new ParallelOptions();
                options.MaxDegreeOfParallelism = coreCount;

                Parallel.For(0, N, options, i =>
                {
                    stock_price[i, 0] = S;
                    stock_price_max[i] = S;
                    stock_price_min[i] = S;

                    for (int j = 1; j < steps; j++)
                    {
                        local_t[i, j] = j * dt;
                        delta[i, j] = option_mc_delta(stock_price[i, j - 1], K, T - local_t[i, j - 1], r, sigma, call_put);
                        stock_price[i, j] = stock_price[i, j - 1] * Math.Exp(nudt + sigmadt * seq[i, j]);
                        cv[i, j] = delta[i, j] * (stock_price[i, j] - stock_price[i, j - 1] * erddt);

                        stock_price_max[i] = Math.Max(stock_price_max[i], stock_price[i, j]);
                        stock_price_min[i] = Math.Min(stock_price_min[i], stock_price[i, j]);
                    }

                    if (call_put == "call")    
                    {
                        CT[i] = Math.Max(0, stock_price_max[i] - K) + beta1 * cv[i, steps - 1];
                    }

                    else                  
                    {
                        CT[i] = Math.Max(0, K - stock_price_min[i]) + beta1 * cv[i, steps - 1];
                    }
                });
            }

            price = CT.Average() * Math.Exp(-r * T);
            standard_error = Math.Sqrt(CT.Select(x => Math.Pow(x - CT.Average(), 2)).Sum() / (N - 1)) / Math.Sqrt(N);

            return (price, standard_error);    
        }

        public static (double delta, double gamma, double vega, double theta, double rho) option_cv_greeks(double S, double K, double T, double r, double sigma, long N, int steps, string call_put, string multi_thread)
        {
            double delta, delta1, delta2, gamma, vega, theta, rho;
            delta = (option_cv(S * 1.05, K, T, r, sigma, N, steps, call_put, multi_thread).option_value - option_cv(S * 0.95, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (2 * S * 0.05);

            delta1 = (option_cv(S * 1.05, K, T, r, sigma, N, steps, call_put, multi_thread).option_value - option_cv(S, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (S * 0.05);
            delta2 = (option_cv(S, K, T, r, sigma, N, steps, call_put, multi_thread).option_value - option_cv(S * 0.95, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (S * 0.05);
            gamma = 2 * (delta1 - delta2) / (2 * S * 0.05);

            vega = (option_cv(S, K, T, r, sigma * 1.05, N, steps, call_put, multi_thread).option_value - option_cv(S, K, T, r, sigma * 0.95, N, steps, call_put, multi_thread).option_value) / (sigma * 1.05 - sigma * 0.95);
            theta = (option_cv(S, K, T * 0.95, r, sigma, N, steps, call_put, multi_thread).option_value - option_cv(S, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (T * 0.05);
            rho = (option_cv(S, K, T, r * 1.05, sigma, N, steps, call_put, multi_thread).option_value - option_cv(S, K, T, r * 0.95, sigma, N, steps, call_put, multi_thread).option_value) / (r * 1.05 - r * 0.95);

            return (delta, gamma, vega, theta, rho);
        }


        public static (double option_value, double standard_error) option_anti(double S, double K, double T, double r, double sigma, long N, int steps, string call_put, string multi_thread)
        {
            double price = 0;
            double standard_error = 0;
            double dt = T / steps;

            double[] option_value = new double[N];     
            double[] option_value_anti = new double[N];
            double[,] stock_price = new double[N, steps];    
            double[,] stock_price_anti = new double[N, steps];
            double[,] seq = new double[N, steps];
            double[] stock_price_max = new double[N];
            double[] stock_price_min = new double[N];
            double[] stock_price_anti_max = new double[N];
            double[] stock_price_anti_min = new double[N];

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < steps; j++)
                {
                    seq[i, j] = Random_Number.Box_Muller_Transform();
                }
            }

            if (multi_thread == "no")
            {
                for (int i = 0; i < N; i++)
                {
                    stock_price[i, 0] = S;
                    stock_price_anti[i, 0] = S;
                    stock_price_max[i] = S;
                    stock_price_min[i] = S;
                    stock_price_anti_max[i] = S;
                    stock_price_anti_min[i] = S;

                    for (int j = 1; j < steps; j++)
                    {
                        stock_price[i, j] = stock_price[i, j - 1] * Math.Exp((r - 0.5 * sigma * sigma) * dt + sigma * Math.Sqrt(dt) * seq[i, j]);
                        stock_price_anti[i, j] = stock_price_anti[i, j - 1] * Math.Exp((r - 0.5 * sigma * sigma) * dt - sigma * Math.Sqrt(dt) * seq[i, j]);

                        stock_price_max[i] = Math.Max(stock_price_max[i], stock_price[i, j]);
                        stock_price_min[i] = Math.Min(stock_price_min[i], stock_price[i, j]);
                        stock_price_anti_max[i] = Math.Max(stock_price_anti_max[i], stock_price_anti[i, j]);
                        stock_price_anti_min[i] = Math.Min(stock_price_anti_min[i], stock_price_anti[i, j]);
                    }

                    if (call_put == "call")    
                    {
                        option_value[i] = Math.Exp(-r * T) * Math.Max(stock_price_max[i] - K, 0);
                        option_value_anti[i] = Math.Exp(-r * T) * (Math.Max(stock_price_max[i] - K, 0) + Math.Max(stock_price_anti_max[i] - K, 0)) / 2;
                    }

                    else                  
                    {
                        option_value[i] = Math.Exp(-r * T) * Math.Max(K - stock_price_min[i], 0);
                        option_value_anti[i] = Math.Exp(-r * T) * (Math.Max(K - stock_price_min[i], 0) + Math.Max(K - stock_price_anti_min[i], 0)) / 2;
                    }
                }
            }
            else
            {
                int coreCount = Environment.ProcessorCount;
                ParallelOptions options = new ParallelOptions();
                options.MaxDegreeOfParallelism = coreCount;

                Parallel.For(0, N, options, i =>
                {
                    stock_price[i, 0] = S;
                    stock_price_anti[i, 0] = S;
                    stock_price_max[i] = S;
                    stock_price_min[i] = S;
                    stock_price_anti_max[i] = S;
                    stock_price_anti_min[i] = S;

                    for (int j = 1; j < steps; j++)
                    {
                        stock_price[i, j] = stock_price[i, j - 1] * Math.Exp((r - 0.5 * sigma * sigma) * dt + sigma * Math.Sqrt(dt) * seq[i, j]);
                        stock_price_anti[i, j] = stock_price_anti[i, j - 1] * Math.Exp((r - 0.5 * sigma * sigma) * dt - sigma * Math.Sqrt(dt) * seq[i, j]);

                        stock_price_max[i] = Math.Max(stock_price_max[i], stock_price[i, j]);
                        stock_price_min[i] = Math.Min(stock_price_min[i], stock_price[i, j]);
                        stock_price_anti_max[i] = Math.Max(stock_price_anti_max[i], stock_price_anti[i, j]);
                        stock_price_anti_min[i] = Math.Min(stock_price_anti_min[i], stock_price_anti[i, j]);
                    }

                    if (call_put == "call")    
                    {
                        option_value[i] = Math.Exp(-r * T) * Math.Max(stock_price_max[i] - K, 0);
                        option_value_anti[i] = Math.Exp(-r * T) * (Math.Max(stock_price_max[i] - K, 0) + Math.Max(stock_price_anti_max[i] - K, 0)) / 2;
                    }

                    else                  
                    {
                        option_value[i] = Math.Exp(-r * T) * Math.Max(K - stock_price_min[i], 0);
                        option_value_anti[i] = Math.Exp(-r * T) * (Math.Max(K - stock_price_min[i], 0) + Math.Max(K - stock_price_anti_min[i], 0)) / 2;
                    }
                });
            }

            price = option_value_anti.Average();
            standard_error = Math.Sqrt(option_value_anti.Select(x => Math.Pow(x - option_value_anti.Average(), 2)).Sum() / (N - 1)) / Math.Sqrt(N);

            return (price, standard_error);
        }

        public static (double delta, double gamma, double vega, double theta, double rho) option_anti_greeks(double S, double K, double T, double r, double sigma, long N, int steps, string call_put, string multi_thread)
        {
            double delta, delta1, delta2, gamma, vega, theta, rho;
            delta = (option_anti(S * 1.05, K, T, r, sigma, N, steps, call_put, multi_thread).option_value - option_anti(S * 0.95, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (2 * S * 0.05);

            delta1 = (option_anti(S * 1.05, K, T, r, sigma, N, steps, call_put, multi_thread).option_value - option_anti(S, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (S * 0.05);
            delta2 = (option_anti(S, K, T, r, sigma, N, steps, call_put, multi_thread).option_value - option_anti(S * 0.95, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (S * 0.05);
            gamma = 2 * (delta1 - delta2) / (2 * S * 0.05);

            vega = (option_anti(S, K, T, r, sigma * 1.05, N, steps, call_put, multi_thread).option_value - option_anti(S, K, T, r, sigma * 0.95, N, steps, call_put, multi_thread).option_value) / (sigma * 1.05 - sigma * 0.95);
            theta = (option_anti(S, K, T * 0.95, r, sigma, N, steps, call_put, multi_thread).option_value - option_anti(S, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (T * 0.05);
            rho = (option_anti(S, K, T, r * 1.05, sigma, N, steps, call_put, multi_thread).option_value - option_anti(S, K, T, r * 0.95, sigma, N, steps, call_put, multi_thread).option_value) / (r * 1.05 - r * 0.95);

            return (delta, gamma, vega, theta, rho);
        }

         public static (double option_value, double standard_error) option_anti_cv(double S, double K, double T, double r, double sigma, long N, int steps, string call_put, string multi_thread)
        {
            double price = 0;
            double standard_error = 0;
            double dt = T / steps;

            double nudt = (r - 0.5 * sigma * sigma) * dt;        
            double sigmadt = sigma * Math.Sqrt(dt);
            double erddt = Math.Exp(r * dt);
            double beta1 = -1;

            double[] option_value = new double[N];     
            double[,] stock_price_1 = new double[N, steps];    
            double[,] stock_price_2 = new double[N, steps];
            double[,] seq = new double[N, steps];
            double[,] cv_1 = new double[N, steps];
            double[,] cv_2 = new double[N, steps];
            double[] CT = new double[N];
            double[,] local_t = new double[N, steps];
            double[,] delta_1 = new double[N, steps];
            double[,] delta_2 = new double[N, steps];
            double[] stock_price_1_max = new double[N];
            double[] stock_price_1_min = new double[N];
            double[] stock_price_2_max = new double[N];
            double[] stock_price_2_min = new double[N];

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < steps; j++)
                {
                    seq[i, j] = Random_Number.Box_Muller_Transform();
                }
            }

            if (multi_thread == "no")
            {
                for (int i = 0; i < N; i++)
                {
                    stock_price_1[i, 0] = S;
                    stock_price_2[i, 0] = S;
                    stock_price_1_max[i] = S;
                    stock_price_1_min[i] = S;
                    stock_price_2_max[i] = S;
                    stock_price_2_min[i] = S;


                    for (int j = 1; j < steps; j++)
                    {
                        local_t[i, j] = j * dt;
                        delta_1[i, j] = option_mc_delta(stock_price_1[i, j - 1], K, T - local_t[i, j - 1], r, sigma, call_put);
                        delta_2[i, j] = option_mc_delta(stock_price_2[i, j - 1], K, T - local_t[i, j - 1], r, sigma, call_put);
                        stock_price_1[i, j] = stock_price_1[i, j - 1] * Math.Exp(nudt + sigmadt * seq[i, j]);
                        stock_price_2[i, j] = stock_price_2[i, j - 1] * Math.Exp(nudt - sigmadt * seq[i, j]);
                        cv_1[i, j] = delta_1[i, j] * (stock_price_1[i, j] - stock_price_1[i, j - 1] * erddt);
                        cv_2[i, j] = delta_2[i, j] * (stock_price_2[i, j] - stock_price_2[i, j - 1] * erddt);

                        stock_price_1_max[i] = Math.Max(stock_price_1_max[i], stock_price_1[i, j]);
                        stock_price_1_min[i] = Math.Min(stock_price_1_min[i], stock_price_1[i, j]);
                        stock_price_2_max[i] = Math.Max(stock_price_2_max[i], stock_price_2[i, j]);
                        stock_price_2_min[i] = Math.Min(stock_price_2_min[i], stock_price_2[i, j]);
                    }

                    if (call_put == "call")    
                    {
                        CT[i] = 0.5 * (Math.Max(0, stock_price_1_max[i] - K) + beta1 * cv_1[i, steps - 1] + Math.Max(0, stock_price_2_max[i] - K) + beta1 * cv_2[i, steps - 1]);
                    }

                    else                  
                    {
                        CT[i] = 0.5 * (Math.Max(0, K - stock_price_1_min[i]) + beta1 * cv_1[i, steps - 1] + Math.Max(0, K - stock_price_2_min[i]) + beta1 * cv_2[i, steps - 1]);
                    }
                }
            }
            else
            {
                int coreCount = Environment.ProcessorCount;
                ParallelOptions options = new ParallelOptions();
                options.MaxDegreeOfParallelism = coreCount;

                Parallel.For(0, N, options, i =>
                {
                    stock_price_1[i, 0] = S;
                    stock_price_2[i, 0] = S;
                    stock_price_1_max[i] = S;
                    stock_price_1_min[i] = S;
                    stock_price_2_max[i] = S;
                    stock_price_2_min[i] = S;

                    for (int j = 1; j < steps; j++)
                    {
                        local_t[i, j] = j * dt;
                        delta_1[i, j] = option_mc_delta(stock_price_1[i, j - 1], K, T - local_t[i, j - 1], r, sigma, call_put);
                        delta_2[i, j] = option_mc_delta(stock_price_2[i, j - 1], K, T - local_t[i, j - 1], r, sigma, call_put);
                        stock_price_1[i, j] = stock_price_1[i, j - 1] * Math.Exp(nudt + sigmadt * seq[i, j]);
                        stock_price_2[i, j] = stock_price_2[i, j - 1] * Math.Exp(nudt - sigmadt * seq[i, j]);
                        cv_1[i, j] = delta_1[i, j] * (stock_price_1[i, j] - stock_price_1[i, j - 1] * erddt);
                        cv_2[i, j] = delta_2[i, j] * (stock_price_2[i, j] - stock_price_2[i, j - 1] * erddt);

                        stock_price_1_max[i] = Math.Max(stock_price_1_max[i], stock_price_1[i, j]);
                        stock_price_1_min[i] = Math.Min(stock_price_1_min[i], stock_price_1[i, j]);
                        stock_price_2_max[i] = Math.Max(stock_price_2_max[i], stock_price_2[i, j]);
                        stock_price_2_min[i] = Math.Min(stock_price_2_min[i], stock_price_2[i, j]);
                    }

                    if (call_put == "call")    
                    {
                        CT[i] = 0.5 * (Math.Max(0, stock_price_1_max[i] - K) + beta1 * cv_1[i, steps - 1] + Math.Max(0, stock_price_2_max[i] - K) + beta1 * cv_2[i, steps - 1]);
                    }

                    else                  
                    {
                        CT[i] = 0.5 * (Math.Max(0, K - stock_price_1_min[i]) + beta1 * cv_1[i, steps - 1] + Math.Max(0, K - stock_price_2_min[i]) + beta1 * cv_2[i, steps - 1]);
                    }
                });
            }

            price = CT.Average() * Math.Exp(-r * T);
            standard_error = Math.Sqrt(CT.Select(x => Math.Pow(x - CT.Average(), 2)).Sum() / (N - 1)) / Math.Sqrt(N);

            return (price, standard_error);
        }

        public static (double delta, double gamma, double vega, double theta, double rho) option_anti_cv_greeks(double S, double K, double T, double r, double sigma, long N, int steps, string call_put, string multi_thread)
        {
            double delta, delta1, delta2, gamma, vega, theta, rho;
            delta = (option_anti_cv(S * 1.05, K, T, r, sigma, N, steps, call_put, multi_thread).option_value - option_anti_cv(S * 0.95, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (2 * S * 0.05);

            delta1 = (option_anti_cv(S * 1.05, K, T, r, sigma, N, steps, call_put, multi_thread).option_value - option_anti_cv(S, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (S * 0.05);
            delta2 = (option_anti_cv(S, K, T, r, sigma, N, steps, call_put, multi_thread).option_value - option_anti_cv(S * 0.95, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (S * 0.05);
            gamma = 2 * (delta1 - delta2) / (2 * S * 0.05);

            vega = (option_anti_cv(S, K, T, r, sigma * 1.05, N, steps, call_put, multi_thread).option_value - option_anti_cv(S, K, T, r, sigma * 0.95, N, steps, call_put, multi_thread).option_value) / (sigma * 1.05 - sigma * 0.95);
            theta = (option_anti_cv(S, K, T * 0.95, r, sigma, N, steps, call_put, multi_thread).option_value - option_anti_cv(S, K, T, r, sigma, N, steps, call_put, multi_thread).option_value) / (T * 0.05);
            rho = (option_anti_cv(S, K, T, r * 1.05, sigma, N, steps, call_put, multi_thread).option_value - option_anti_cv(S, K, T, r * 0.95, sigma, N, steps, call_put, multi_thread).option_value) / (r * 1.05 - r * 0.95);

            return (delta, gamma, vega, theta, rho);
        }
    }
}
