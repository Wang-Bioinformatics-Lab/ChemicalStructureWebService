

def rdkit_handle_error(func):
    """Run function as normal - catch TypeError from RDKit caused by C++ error
    Server will still revert to 500 status when other error occur.
    """
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except (TypeError, NameError) as e:
            print("RDKit unable to process mol")
            #print(e, file=sys.stderr, flush=True)
            return {"message": "structure cant be identified"}, 400
    # Have to rename wrapper or else Flask tries to register all uses as "wrapper"
    # causing AssertionError
    wrapper.__name__ = func.__name__
    return wrapper