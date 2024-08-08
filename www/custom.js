/*
* Observer for canvasXpress output to be displayed or
* all plots are plots are unavailable
*/
function waitForElement(selector, matrix_size) {
    return new Promise(resolve => {
        if (document.querySelector(selector)) {
            return resolve(document.querySelector(selector));
        }

        const observer = new MutationObserver(mutations => {
            no_of_not_available_plots = $('#matrixPlotOutput').find('.text-danger').length;

            if ((no_of_not_available_plots == matrix_size) ||
                (document.querySelector(selector))) {
                    resolve(document.querySelector(selector));
                    observer.disconnect();
            }
        });

        observer.observe(document.body, {
            childList: true,
            subtree: true
        });
    });
}


/*
* Dismiss matrix plot render modal
*/
function dismissMatrixBusyModal(matrix_size) {
    waitForElement('#matrixPlotOutput canvas', matrix_size).then((element) => {
        Shiny.setInputValue('matrixRenderIsComplete', 'TRUE', {priority: 'event'});
   })

}


/* hide datatable filter row when table is empty or all filter inputs are disabled*/
function hideSearchRow() {
   $(document).on("xhr.dt", function(e) {
       if ($(".dataTables_scrollHead").find("input.form-control[disabled]").length == $(".dataTables_scrollHead").find("input.form-control").length) {
           $("div.dataTables_wrapper > div.dataTables_scroll > div.dataTables_scrollHead > div > table > thead > tr:nth-child(2)").hide();
       }
   });
}


// make actionButton click independent of row selection in table
function stopClickPropagation(event) {
    event.stopPropagation();
}
