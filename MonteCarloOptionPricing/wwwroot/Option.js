document.getElementById("button1").addEventListener("click", function () { MonteCarlo_sim() });


function read_ui_inputs() {
    var inputs = {};
    inputs.Call_put = document.getElementById("Call_put").value;
    inputs.S = document.getElementById("S").value;
    inputs.K = document.getElementById("K").value;
    inputs.T = document.getElementById("T").value;
    inputs.r = document.getElementById("r").value;
    inputs.sigma = document.getElementById("sigma").value;
    inputs.N = document.getElementById("N").value;
    inputs.steps = document.getElementById("steps").value;
    inputs.Anti = document.getElementById("Anti").value;
    inputs.CV = document.getElementById("CV").value;
    inputs.Multi_thread = document.getElementById("Multi_thread").value;
    inputs.Option_type = document.getElementById("Option_type").value;
    inputs.pay_amount = document.getElementById("pay_amount").value;
    inputs.Barrier_type = document.getElementById("Barrier_type").value;
    inputs.barrier_level = document.getElementById("barrier_level").value;
    return inputs;
}



function MonteCarlo_sim() {
    console.log(1)
    var inputs = read_ui_inputs();
    console.log(2)

    let request = new XMLHttpRequest();
    console.log(3)
    request.open("POST", "http://localhost:5094/Main", true);
    console.log(4)
    request.setRequestHeader("Access-Control-Allow-Origin", "*");
    console.log(5)
    request.setRequestHeader("Accept", "application/json");
    request.setRequestHeader("Content-Type", "application/json;charset=UTF-8");


    request.onreadystatechange = () => {
        console.log("error!");
        if (request.readyState == XMLHttpRequest.DONE && request.status == 200) {
            var results = JSON.parse(request.response);
            document.getElementById("option_value").innerHTML = results.option_value;
            document.getElementById("standard_error").innerHTML = results.standard_error;
            document.getElementById("delta").innerHTML = results.delta;
            document.getElementById("gamma").innerHTML = results.gamma;
            document.getElementById("vega").innerHTML = results.vega;
            document.getElementById("theta").innerHTML = results.theta;
            document.getElementById("rho").innerHTML = results.rho;
        }
        else if (request.readyState == XMLHttpRequest.DONE && request.status == 400) {
            console.log("error!");
        }
        else {
            console.log("error!!!");
        }
    };

    request.send(JSON.stringify(inputs));

}
