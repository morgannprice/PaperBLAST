
// to be used as an onClick method
// (Figures out what the target URL is from the item parameter)
function logger(item, label) {
    var URL = item.href;
    var xhttp = new XMLHttpRequest();
    // use a random code to prevent caching
    var request = "../cgi-bin/noop.cgi?label=" + label + "&URL=" + encodeURI(URL) + "&code=" + Math.random();
    xhttp.open("GET", request);
    xhttp.send();
    return(true);
}
