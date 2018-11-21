
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

// for expanding rows in a table -- assumes that the specified row contains the expand link or button
//	so it should be hidden, and that the very 1st row of the table is always shown
function tblexpander(obj) {
    var tr = obj;
    while(tr && tr.tagName != "TR") {
      tr = tr.parentNode;
    }
    var trlist = tr.parentNode.children;
    var i;
    // assume the 1st row already has the desired display value
    for (i = 1; i < trlist.length; i++) {
        trlist[i].style.display = trlist[0].style.display;
    }
    tr.style.display = "none";
    return(true);
}
