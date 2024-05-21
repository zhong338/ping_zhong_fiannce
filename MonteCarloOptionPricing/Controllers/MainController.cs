// This file is in collaboration with my classmate Zongxing Li. 
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.SignalR;
using MonteCarloOptionPricing.Options;

namespace MonteCarloOptionPricing.Controllers
{
    [ApiController]
    [Route("[controller]")]
    public class MainController : ControllerBase
    {
        public class Simulation_inputs
        {
            public double S { get; set; }
            public double K { get; set; }
            public double T { get; set; }
            public double r { get; set; }
            public double sigma { get; set; }
            public long N { get; set; }
            public int steps { get; set; }
            public string? Call_put { get; set; }
            public string? Option_type { get; set; }
            public string? Anti { get; set; }
            public string? CV { get; set; }
            public string? Multi_thread { get; set; }
            public double pay_amount { get; set; }
            public double barrier_level { get; set; }
            public string? Barrier_type { get; set; }
        }

        public class Simulation_results
        {
            public double option_value { get; set; }
            public double standard_error { get; set; }
            public double delta { get; set; }
            public double gamma { get; set; }
            public double vega { get; set; }
            public double theta { get; set; }
            public double rho { get; set; }
        }


        [HttpPost]
        public ActionResult<Simulation_results> Post([FromBody] Simulation_inputs inputs)
        {
            double S = inputs.S;
            double K = inputs.K;
            double T = inputs.T;
            double r = inputs.r;
            double sigma = inputs.sigma;
            long N = inputs.N;
            int steps = inputs.steps;
            string call_put = inputs.Call_put!;
            string option_type = inputs.Option_type!;
            string anti = inputs.Anti!;
            string CV = inputs.CV!;
            string multi_thread = inputs.Multi_thread!;
            double pay_amount = inputs.pay_amount;
            double barrier_level = inputs.barrier_level;
            string barrier_type = inputs.Barrier_type!;

            double option_value = 0;
            double standard_error = 0;
            double delta = 0;
            double gamma = 0;
            double vega = 0;
            double theta = 0;
            double rho = 0;

            switch (option_type)
            {
                case "European":      

                    if (anti == "no")
                    {
                        if (CV == "no")
                        {
                            (option_value, standard_error) = EuropeanOption.option_origin(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                            (delta, gamma, vega, theta, rho) = EuropeanOption.option_origin_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                        }
                        else
                        {
                            (option_value, standard_error) = EuropeanOption.option_cv(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                            (delta, gamma, vega, theta, rho) = EuropeanOption.option_cv_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                        }
                    }
                    else
                    {
                        if (CV == "no")
                        {
                            (option_value, standard_error) = EuropeanOption.option_anti(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                            (delta, gamma, vega, theta, rho) = EuropeanOption.option_anti_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                        }
                        else
                        {
                            (option_value, standard_error) = EuropeanOption.option_anti_cv(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                            (delta, gamma, vega, theta, rho) = EuropeanOption.option_anti_cv_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                        }
                    }
                    break;


                case "Asian":      

                    if (anti == "no")
                    {
                        if (CV == "no")
                        {
                            (option_value, standard_error) = AsianOption.option_origin(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                            (delta, gamma, vega, theta, rho) = AsianOption.option_origin_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                        }
                        else
                        {
                            (option_value, standard_error) = AsianOption.option_cv(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                            (delta, gamma, vega, theta, rho) = AsianOption.option_cv_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                        }
                    }
                    else
                    {
                        if (CV == "no")
                        {
                            (option_value, standard_error) = AsianOption.option_anti(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                            (delta, gamma, vega, theta, rho) = AsianOption.option_anti_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                        }
                        else
                        {
                            (option_value, standard_error) = AsianOption.option_anti_cv(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                            (delta, gamma, vega, theta, rho) = AsianOption.option_anti_cv_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                        }
                    }
                    break;


                case "Digital":

                    Console.Write("Pay Amount For Digital Option: ");
                    pay_amount = Convert.ToDouble(Console.ReadLine());

                    if (anti == "no")
                    {
                        if (CV == "no")
                        {
                            (option_value, standard_error) = DigitalOption.option_origin(S, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount);
                            (delta, gamma, vega, theta, rho) = DigitalOption.option_origin_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount);
                        }
                        else
                        {
                            Console.WriteLine("");
                            Console.WriteLine("Control Variate is not available for Digital option.");     
                        }
                    }
                    else
                    {
                        if (CV == "no")
                        {
                            (option_value, standard_error) = DigitalOption.option_anti(S, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount);
                            (delta, gamma, vega, theta, rho) = DigitalOption.option_anti_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread, pay_amount);
                        }
                        else
                        {
                            Console.WriteLine("");
                            Console.WriteLine("Control Variate is not available for Digital option.");
                        }
                    }

                    break;

                case "Barrier":

                    if (anti == "no")
                    {
                        if (CV == "no")
                        {
                            (option_value, standard_error) = BarrierOption.option_origin(S, K, T, r, sigma, N, steps, call_put, multi_thread, barrier_type, barrier_level);
                            (delta, gamma, vega, theta, rho) = BarrierOption.option_origin_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread, barrier_type, barrier_level);
                        }
                        else
                        {
                            (option_value, standard_error) = BarrierOption.option_cv(S, K, T, r, sigma, N, steps, call_put, multi_thread, barrier_type, barrier_level);
                            (delta, gamma, vega, theta, rho) = BarrierOption.option_cv_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread, barrier_type, barrier_level);
                        }
                    }
                    else
                    {
                        if (CV == "no")
                        {
                            (option_value, standard_error) = BarrierOption.option_anti(S, K, T, r, sigma, N, steps, call_put, multi_thread, barrier_type, barrier_level);
                            (delta, gamma, vega, theta, rho) = BarrierOption.option_anti_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread, barrier_type, barrier_level);
                        }
                        else
                        {
                            (option_value, standard_error) = BarrierOption.anti_cv(S, K, T, r, sigma, N, steps, call_put, multi_thread, barrier_type, barrier_level);
                            (delta, gamma, vega, theta, rho) = BarrierOption.anti_cv_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread, barrier_type, barrier_level);
                        }
                    }
                    break;


                case "Lookback":   
                    if (anti == "no")
                    {
                        if (CV == "no")
                        {
                            (option_value, standard_error) = LookbackOption.option_origin(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                            (delta, gamma, vega, theta, rho) = LookbackOption.option_origin_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                        }
                        else
                        {
                            (option_value, standard_error) = LookbackOption.option_cv(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                            (delta, gamma, vega, theta, rho) = LookbackOption.option_cv_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                        }
                    }
                    else
                    {
                        if (CV == "no")
                        {
                            (option_value, standard_error) = LookbackOption.option_anti(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                            (delta, gamma, vega, theta, rho) = LookbackOption.option_anti_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                        }
                        else
                        {
                            (option_value, standard_error) = LookbackOption.option_anti_cv(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                            (delta, gamma, vega, theta, rho) = LookbackOption.option_anti_cv_greeks(S, K, T, r, sigma, N, steps, call_put, multi_thread);
                        }
                    }
                    break;


                case "Range":       
                    if (anti == "no")
                    {
                        if (CV == "no")
                        {
                            (option_value, standard_error) = RangeOption.option_origin(S, T, r, sigma, N, steps, multi_thread);
                            (delta, gamma, vega, theta, rho) = RangeOption.option_origin_greeks(S, T, r, sigma, N, steps, multi_thread);
                        }
                        else
                        {
                            Console.WriteLine("");
                            Console.WriteLine("Control Variate is not available for Range option.");     
                        }
                    }
                    else
                    {
                        if (CV == "no")
                        {
                            (option_value, standard_error) = RangeOption.option_anti(S, T, r, sigma, N, steps, multi_thread);
                            (delta, gamma, vega, theta, rho) = RangeOption.option_anti_greeks(S, T, r, sigma, N, steps, multi_thread);
                        }
                        else
                        {
                            Console.WriteLine("");
                            Console.WriteLine("Control Variate is not available for Range option.");
                        }
                    }
                    break;
            }

            var results = new Simulation_results();
            results.option_value = option_value;
            results.standard_error = standard_error;

            results.delta = delta;
            results.gamma = gamma;
            results.vega = vega;
            results.theta = theta;
            results.rho = rho;

            return (Ok(results));
        }
    }
}
