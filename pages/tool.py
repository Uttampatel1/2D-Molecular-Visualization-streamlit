import streamlit as st

def main():
    st.title("MolView Embedding Example")

    # Embedding the MolView link using HTML
    molview_iframe = '<iframe style="width: 500px; height: 300px;" frameborder="0" src="https://embed.molview.org"></iframe>'
    st.markdown(molview_iframe, unsafe_allow_html=True)

if __name__ == "__main__":
    main()
