// based on
// https://www.w3schools.com/howto/howto_js_autocomplete.asp

function autocomplete(inp) {
    var arr = [];
    function closeAllLists(elmnt) {
        var x = document.getElementsByClassName("autocomplete-items");
        for (var i = 0; i < x.length; i++) {
            if (elmnt != x[i] && elmnt != inp) {
                x[i].parentNode.removeChild(x[i]); 
            }
        }
        arr = [];
        console.log("emptied");
    }

    var currentFocus;
    inp.addEventListener("input", function(e) {
            var obj = this; // current object
            var val = this.value; // current text
            closeAllLists();
            if (!val) { return false; }
            var URL = "uniprot_proteome_expand.cgi?query=" + encodeURIComponent(val);
            console.log(URL);
            fetch(URL).then(function(response) {
                    if(!response.ok)  {
                        console.log("Cannot contact server");
                    } else {
                        return response.text();
                    }
                }).then(function(text) {
                        if (text) {
                            arr = [];
                            var lines = text.split("\n");
                            for (i = 0; i < lines.length; i++) {
                                var pieces = lines[i].split("\t");
                                if(pieces.length > 0) {
                                    arr.push(pieces[0]);
                                }
                            }
                            console.log("filled to " + arr.length);
                            currentFocus = -1;
                            a = document.createElement("DIV");
                            a.setAttribute("id", obj.id + "autocomplete-list");
                            a.setAttribute("class","autocomplete-items");
                            obj.parentNode.appendChild(a);
                            for (i = 0; i < arr.length; i++) {
                                if(arr[i].substr(0,val.length).toUpperCase() == val.toUpperCase()) {
                                    b = document.createElement("DIV");
                                    b.innerHTML = "<strong>" + arr[i].substr(0,val.length) + "</strong>" + arr[i].substr(val.length) +
                                        "<input type='hidden' value='" + arr[i] + "'>";
                                    b.addEventListener("click", function(e) {
                                            inp.value = this.getElementsByTagName("input")[0].value;
                                            closeAllLists();
                                        });
                                    a.appendChild(b);
                                }
                            }
                        }
                    });
        });
    inp.addEventListener("keydown", function(e) {
            var x = document.getElementById(this.id + "autocomplete-list");
            if (x) x = x.getElementsByTagName("div");
            if (e.keyCode==40) { // down key
                currentFocus++; addActive(x);
            } else if (e.keyCode==38) { // up key
                currentFocus--; addActive(x);
            } else if (e.keyCode==13) { // enter key
                e.preventDefault();
                if (currentFocus > -1) {
                    if(x) x[currentFocus].click();
                }
            }
        });

    function addActive(x) {
        if(!x) return false;
        removeActive(x);
        if (currentFocus >= x.length) currentFocus = 0;
        if (currentFocus < 0) currentFocus = (x.length - 1);
        x[currentFocus].classList.add("autocomplete-active");
    }
   function removeActive(x) {
       for (var i = 0; i < x.length; i++) {
           x[i].classList.remove("autocomplete-active");
       }
   }
   document.addEventListener("click", function (e) { closeAllLists(e.target); });
}


