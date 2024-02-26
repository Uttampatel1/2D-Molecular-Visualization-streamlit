import streamlit as st
import numpy as np

def pic50_to_ic50(pic50, unit):
    if unit == "millimeter":
        ic50 = np.power(10, -pic50)
    elif unit == "micrometer":
        ic50 = np.power(10, -6-pic50)
    elif unit == "nanometer":
        ic50 = np.power(10, -9-pic50)
    return ic50

def ic50_to_pic50(ic50, unit):
    pic50 = None
    if unit == "millimeter":
        pic50 = -np.log10(ic50)
    elif unit == "micrometer":
        pic50 = -np.log10(ic50 * 1e6)
    elif unit == "nanometer":
        pic50 = -np.log10(ic50 * 1e9)
    return pic50

def main():
    st.title("pIC50 to IC50 Converter")

    conversion_type = st.radio("Select conversion type:", ("pIC50 to IC50", "IC50 to pIC50"))

    if conversion_type == "pIC50 to IC50":
        pic50 = st.number_input("Enter pIC50 value:", min_value=0.0, step=0.1)
        unit = st.selectbox("Select unit of concentration:", ("millimeter", "micrometer", "nanometer"))
        if st.button("Convert"):
            ic50 = pic50_to_ic50(pic50, unit)
            st.success(f"IC50 value: {ic50:.6f} {unit}")

    elif conversion_type == "IC50 to pIC50":
        ic50 = st.number_input("Enter IC50 value:", min_value=0.0, step=0.000001)
        unit = st.selectbox("Select unit of concentration:", ("millimeter", "micrometer", "nanometer"))
        if st.button("Convert"):
            pic50 = ic50_to_pic50(ic50, unit)
            st.success(f"pIC50 value: {pic50:.2f}")

if __name__ == "__main__":
    main()
