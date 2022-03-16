
function leafSearch() {
    var nMatch = 0;
    var query = document.getElementById("query").value;
    query = query.replace(/^ +/, "");
    query = query.replace(/ +$/, "");
    if (query == "") return false;
    var queryU = query.toUpperCase();
    var labels = document.getElementById("gLabels");
    if (labels) {
        // find matching label at right
        const list = labels.childNodes;
        for(let i = 0; i < list.length; i++) {
            let o = list[i];
            let title = o.firstChild;
            if (title && matchQuery(queryU, title.textContent)) {
                // setting the style for the whole object does not work in chrome
                // instead, set style for the children after the  title
                // (if there is a link, there is the A tag; otherwise, two tspans)
                o.children[1].style.fill="blue";
                if(o.children.length > 2) o.children[2].style.fill="blue";
                nMatch++;
            }
        }
    } else {
        // no labels, find the matching node and turn the text on
        var gs = document.getElementsByTagName("g");
        for (i = 0; i < gs.length; i++) {
            var list = gs[i].children;
            if (list.length == 2 && list[0].tagName == "circle" && list[1].tagName=="text") {
                var textObject = list[1];
                if (matchQuery(queryU, textObject.textContent)) {
                    textObject.style.display = "inline";
                    textObject.style.fill = "blue";
                    nMatch++;
                }
            }
        }
    }
    document.getElementById("searchStatement").textContent = nMatch + " matched "
        + query + " (including partial matches to words)";
    return false;
}

// match partial words at beginning or after space
function matchQuery(queryU, textContent) {
    var content = textContent.toUpperCase();
    return content.substring(0, queryU.length) == queryU
        || content.indexOf(" " + queryU) >= 0;
}

function leafClear() {
    var labels = document.getElementById("gLabels");
    if (labels) {
        const list = labels.childNodes;
        for(let i = 0; i < list.length; i++) {
            let o = list[i];
            let title = list[i].firstChild;
            if (title) {
                o.children[1].style.removeProperty("fill");
                if(o.children.length > 2) o.children[2].style.removeProperty("fill");
            }
        }
    }
    // and clear popup labels
    var gs = document.getElementsByTagName("g");
    for (i = 0; i < gs.length; i++) {
        var list = gs[i].children;
        if (list.length == 2 && list[0].tagName == "circle" && list[1].tagName=="text") {
            var textObject = list[1];
            textObject.style.display = "none";
        }
    }
    document.getElementById("query").value = "";
    document.getElementById("searchStatement").textContent = "";
    return false;
}

function leafClick(o) {
  var textObject = o.parentNode.children[1];
  textObject.style.display = "inline";
  textObject.style.fill = "#AA0000";
}
