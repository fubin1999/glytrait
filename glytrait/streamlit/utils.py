import streamlit as st


def cache_wrapper(func):
    @st.cache_data(ttl=3660)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


@st.cache_data(ttl=600)
def get_description(trait: str):
    """Get the defination of a derived trait."""
    for f in st.session_state.formulas:
        if f.name == trait:
            return f.description
