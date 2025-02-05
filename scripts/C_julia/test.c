//gcc -o julia_test test.c -I/Applications/Julia-1.10.app/Contents/Resources/julia/include/julia -L/Applications/Julia-1.10.app/Contents/Resources/julia/lib -ljulia
#include <stdio.h>
#include <julia.h>

int main() {
    // Initialize the Julia runtime
    jl_init();

    // Call a simple Julia expression to print "Hello from Julia!"
    jl_eval_string("println(\"Hello from Julia!\")");

    // Call a Julia expression to create a variable 'x'
    jl_eval_string("x = 42");
    jl_value_t *x = jl_get_global(jl_main_module, jl_symbol("x"));

    // Call the Julia 'show' function with 'x' as an argument
    jl_value_t *args[1] = {x}; // Array holding the argument
    jl_call1(jl_get_function(jl_main_module, "show"), args[0]); // Call 'show' with x as the argument

    // Clean up and shut down the Julia runtime
    jl_atexit_hook(0);
    return 0;
}