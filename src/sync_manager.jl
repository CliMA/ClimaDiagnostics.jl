struct NetCDFLockManager
    """A mapping from open NetCDF files to the corresponding lock"""
    locks::Dict{NCDatasets.NCDataset, ReentrantLock}

    """A channel holding all the request for syncing"""
    sync_requests::Channel{Tuple{NCDatasets.NCDataset, ReentrantLock}}

    """A task waiting for sync requests"""
    sync_task::Task
end

function NetCDFLockManager()
    locks = Dict{NCDatasets.NCDataset, ReentrantLock}()
    sync_requests = Channel{Tuple{NCDatasets.NCDataset, ReentrantLock}}()

    # This automatically close when the channel is closed
    # TODO: Look into if isopen is the right thing here?
    sync_task = Threads.@spawn while isopen(sync_requests)
        nc, lk = take!(sync_requests)
        @info "Syncing data on IO thread for $(keys(nc))"
        lock(lk)
        try
            NCDatasets.sync(nc)
        finally
            unlock(lk)
        end
        @info "Finished syncing data on IO thread for $(keys(nc))"
    end

    return NetCDFLockManager(locks, sync_requests, sync_task)
end

function get_lock(lock_manager::NetCDFLockManager, nc::NCDatasets.NCDataset)
    lk = get!(lock_manager.locks, nc) do
        ReentrantLock()
    end
    return lk
end

function add_sync_request(
    lock_manager::NetCDFLockManager,
    nc::NCDatasets.NCDataset,
)
    lk = get_lock(lock_manager, nc)
    put!(lock_manager.sync_requests, (nc, lk))
    return nothing
end

function Base.close(lock_manager::NetCDFLockManager)
    close(lock_manager.sync_requests)
    return nothing
end
