
shinyjs.moveButton = function(params) {
    $('#' + params[0]).append($('#' + params[1]));
}


shinyjs.removeInputElems = function() {
    $('div[style*="display: flex;"]').remove()
}