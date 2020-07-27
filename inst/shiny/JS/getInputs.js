
shinyjs.getInputElems = function(params){
    var inputelems = $('select[id^=' + params + ']');
    var inputids = [];
    l = inputelems.length;
    for (i = 0; i < l; i++) {
        inputids.push(inputelems[i].id);
    }
    Shiny.setInputValue(params + 'inputElems',inputids);
};