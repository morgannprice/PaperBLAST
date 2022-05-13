
function getTitleObjects() {
    var labels = document.getElementById("gLabels");
    var alnRows = document.getElementById("alnRows");
    var titles = [];
    if (labels) {
        // find matching label at right in choose-sites mode
        for(const o of labels.children) {
            titles.push(o);
        }
    } else if (alnRows) {
        // find matching label at right in auto-sites mode
        const list = alnRows.children;
        for (const row of alnRows.children) {
            // the first child is the rect for highlighting the row; the second is the text object
            titles.push(row.children[1]);
        }
    } else {
        // no labels, find the text node after each circle
        const gs = document.getElementsByTagName("g");
        for (var g of gs) {
            var list = g.children;
            // the hidden text is after the circle
            if (list.length == 2 && list[0].tagName == "circle" && list[1].tagName=="text") {
                titles.push(list[1]);
            }
        }
    }
    return titles;
}

function leafSearch() {
    var labels = document.getElementById("gLabels");
    var alnRows = document.getElementById("alnRows");
    var query = document.getElementById("query").value;
    query = query.replace(/^ +/, "");
    query = query.replace(/ +$/, "");
    if (query == "") return false;
    var queryU = query.toUpperCase();
    var nMatch = 0;
    for (var title of getTitleObjects()) {
        if (matchQuery(queryU, title.textContent)) {
            if (!labels && !alnRows) {
                title.style.display = "inline";
            }
            // This looks good, but does not work reliably in chrome, not sure why not
            title.style.fill = "blue";
            // Works reliably, a bit funny looking
            title.style.stroke = "blue";
            title.style.strokeWidth = "0.5px";
            nMatch++;
        }
    }
    document.getElementById("searchStatement").textContent = nMatch + " matched "
        + query + " (including partial matches to words)";
    return false;
}

// match partial words at beginning or after space or ( or | or =
// (which is sometimes used as a separator in definition lines)
function matchQuery(queryU, textContent) {
    if (!textContent) { return false; }
    var content = textContent.toUpperCase();
    return content.substring(0, queryU.length) == queryU
        || content.indexOf(" " + queryU) >= 0
        || content.indexOf("|" + queryU) >= 0
        || content.indexOf("=" + queryU) >= 0
        || content.indexOf("(" + queryU) >= 0
        || content.indexOf("_" + queryU) >= 0
        || content.indexOf("[" + queryU) >= 0;
}

function leafClear() {
    var labels = document.getElementById("gLabels");
    var alnRows = document.getElementById("alnRows");
    var titles = getTitleObjects();
    for (var t of getTitleObjects()) {
        t.style.stroke = null;
        t.style.fill = null;
        if (!labels && !alnRows) {
            t.style.display = "none";
        }
    }
    document.getElementById("query").value = "";
    document.getElementById("searchStatement").textContent = "";
    return false;
}

function leafClick(o) {
  var textObject = o.parentNode.children[1];
  textObject.style.display = "inline";
  textObject.style.fill = "blue";
  textObject.style.stroke = "blue";
  textObject.style.strokeWidth = 0.5;
}
