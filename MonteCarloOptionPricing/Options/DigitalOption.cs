using System;
using System.Diagnostics;

// This file is in collaboration with my classmate Zongxing Li. 

namespace MonteCarloOptionPricing.Options
{
    public class DigitalOption
    {

        public static (double option_value, double standard_error) option_origin(double S, double K, double T, double r, double sigma, long N, int steps, string call_put, string multi_thread, double pay_amount)
        {
            double price = 0;
            double standard_error = 0;
            double dt = T / steps;

            double[] option_value = new double[N];     
            double[,] stock_price = new double[N, steps];    
            double[,] seq = new double[N, steps];

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

                    for (int j = 1; j < steps; j++)
                    {
                        stock_price[i, j] = stock_price[i, j - 1] * Math.Exp((r - 0.5 * sigma * sigma) * dt + sigma * Math.Sqrt(dt) * seq[i, j]);
                    }
                    if (call_put == "call")    
                    {
                        option_value[i] = Math.Exp(-r * T) * (stock_price[i, steps - 1] - K > 0 ? pay_amount : 0);  // Conditional operator used to determine payoff
                    }
                    else                  
                    {
                        option_value[i] = Math.Exp(-r * T) * (stock_price[i, steps - 1] - K < 0 ? pay_amount : 0);
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

                    for (int j = 1; j < steps; j++)
                    {
                        stock_price[i, j] = stock_price[i, j - 1] * Math.Exp((r - 0.5 * sigma * sigma) * dt + sigma * Math.Sqrt(dt) * seq[i, j]);
                    }
                    if (call_put == "call")    
                    {
                        option_value[i] = Math.Exp(-r * T) * (stock_price[i, steps - 1] - K > 0 ? pay_amount : 0);
                    }
                    else                  
                    {
                        option_value[i] = Math.Exp(-r * T) * (stock_price[i, steps - 1] - K < 0 ? pay_amount : 0);
                    }
                });
            }

            price = option_value.Average();
            standard_error = Math.Sqrt(option_value.Select(x => Math.Pow(x - option_value.Average(), 2)).Sum() / (N - 1)) / Math.Sqrt(N);

            return (price, standard_error);    
        }

        public static (double delta, double gamma, double vega, double theta, double rho) option_origin_greeks(double S, double K, double T, double r, double sigma, long N, int steps, string call_put, string multi_thread, double pay_amount)
        {
            double delta, delta1, delta2, gamma, vega, theta, rho;
            delta = (option_origin(S * 1.05, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value - option_origin(S * 0.95, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value) / (2 * S * 0.05);

            delta1 = (option_origin(S * 1.05, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value - option_origin(S, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value) / (S * 0.05);
            delta2 = (option_origin(S, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value - option_origin(S * 0.95, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value) / (S * 0.05);
            gamma = 2 * (delta1 - delta2) / (2 * S * 0.05);

            vega = (option_origin(S, K, T, r, sigma * 1.05, N, steps, call_put, multi_thread, pay_amount).option_value - option_origin(S, K, T, r, sigma * 0.95, N, steps, call_put, multi_thread, pay_amount).option_value) / (sigma * 1.05 - sigma * 0.95);
            theta = (option_origin(S, K, T * 0.95, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value - option_origin(S, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value) / (T * 0.05);
            rho = (option_origin(S, K, T, r * 1.05, sigma, N, steps, call_put, multi_thread, pay_amount).option_value - option_origin(S, K, T, r * 0.95, sigma, N, steps, call_put, multi_thread, pay_amount).option_value) / (r * 1.05 - r * 0.95);

            return (delta, gamma, vega, theta, rho);
        }

        

       
        public static (double option_value, double standard_error) option_anti(double S, double K, double T, double r, double sigma, long N, int steps, string call_put, string multi_thread, double pay_amount)
        {
            double price = 0;
            double standard_error = 0;
            double dt = T / steps;

            double[] option_value = new double[N];     
            double[] option_value_anti = new double[N];
            double[,] stock_price = new double[N, steps];    
            double[,] stock_price_anti = new double[N, steps];
            double[,] seq = new double[N, steps];

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

                    for (int j = 1; j < steps; j++)
                    {
                        stock_price[i, j] = stock_price[i, j - 1] * Math.Exp((r - 0.5 * sigma * sigma) * dt + sigma * Math.Sqrt(dt) * seq[i, j]);
                        stock_price_anti[i, j] = stock_price_anti[i, j - 1] * Math.Exp((r - 0.5 * sigma * sigma) * dt - sigma * Math.Sqrt(dt) * seq[i, j]);
                    }

                    if (call_put == "call")    
                    {
                        option_value[i] = Math.Exp(-r * T) * (stock_price[i, steps - 1] - K > 0 ? pay_amount : 0);
                        option_value_anti[i] = Math.Exp(-r * T) * ((stock_price[i, steps - 1] - K > 0 ? pay_amount : 0) + (stock_price_anti[i, steps - 1] - K > 0 ? pay_amount : 0)) / 2;
                    }

                    else                 
                    {
                        option_value[i] = Math.Exp(-r * T) * (stock_price[i, steps - 1] - K < 0 ? pay_amount : 0);
                        option_value_anti[i] = Math.Exp(-r * T) * ((stock_price[i, steps - 1] - K < 0 ? pay_amount : 0) + (stock_price_anti[i, steps - 1] - K < 0 ? pay_amount : 0)) / 2;
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

                    for (int j = 1; j < steps; j++)
                    {
                        stock_price[i, j] = stock_price[i, j - 1] * Math.Exp((r - 0.5 * sigma * sigma) * dt + sigma * Math.Sqrt(dt) * seq[i, j]);
                        stock_price_anti[i, j] = stock_price_anti[i, j - 1] * Math.Exp((r - 0.5 * sigma * sigma) * dt - sigma * Math.Sqrt(dt) * seq[i, j]);
                    }

                    if (call_put == "call")    
                    {
                        option_value[i] = Math.Exp(-r * T) * (stock_price[i, steps - 1] - K > 0 ? pay_amount : 0);
                        option_value_anti[i] = Math.Exp(-r * T) * ((stock_price[i, steps - 1] - K > 0 ? pay_amount : 0) + (stock_price_anti[i, steps - 1] - K > 0 ? pay_amount : 0)) / 2;
                    }

                    else                  
                    {
                        option_value[i] = Math.Exp(-r * T) * (stock_price[i, steps - 1] - K < 0 ? pay_amount : 0);
                        option_value_anti[i] = Math.Exp(-r * T) * ((stock_price[i, steps - 1] - K < 0 ? pay_amount : 0) + (stock_price_anti[i, steps - 1] - K < 0 ? pay_amount : 0)) / 2;
                    }
                });
            }

            price = option_value_anti.Average();
            standard_error = Math.Sqrt(option_value_anti.Select(x => Math.Pow(x - option_value_anti.Average(), 2)).Sum() / (N - 1)) / Math.Sqrt(N);

            return (price, standard_error);
        }

        public static (double delta, double gamma, double vega, double theta, double rho) option_anti_greeks(double S, double K, double T, double r, double sigma, long N, int steps, string call_put, string multi_thread, double pay_amount)
        {
            double delta, delta1, delta2, gamma, vega, theta, rho;
            delta = (option_anti(S * 1.05, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value - option_anti(S * 0.95, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value) / (2 * S * 0.05);

            delta1 = (option_anti(S * 1.05, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value - option_anti(S, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value) / (S * 0.05);
            delta2 = (option_anti(S, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value - option_anti(S * 0.95, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value) / (S * 0.05);
            gamma = 2 * (delta1 - delta2) / (2 * S * 0.05);

            vega = (option_anti(S, K, T, r, sigma * 1.05, N, steps, call_put, multi_thread, pay_amount).option_value - option_anti(S, K, T, r, sigma * 0.95, N, steps, call_put, multi_thread, pay_amount).option_value) / (sigma * 1.05 - sigma * 0.95);
            theta = (option_anti(S, K, T * 0.95, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value - option_anti(S, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount).option_value) / (T * 0.05);
            rho = (option_anti(S, K, T, r * 1.05, sigma, N, steps, call_put, multi_thread, pay_amount).option_value - option_anti(S, K, T, r * 0.95, sigma, N, steps, call_put, multi_thread, pay_amount).option_value) / (r * 1.05 - r * 0.95);

            return (delta, gamma, vega, theta, rho);
        }
    }
}
