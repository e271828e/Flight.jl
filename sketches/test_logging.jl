using Logging

old_logger = global_logger(ConsoleLogger(Logging.Info))
disable_logging(Logging.Warn) #nothing below Warn neither prints nor allocates
# disable_logging(Logging.Debug) #every message allocates, even if not printed

function test_logging()

    @debug "Hello, debug"
    @info "Hello, info"
    @warn "Hello, warn"
    @error "Hello, error"
end