import streamlit as st
import math

def pic50_to_ic50(pic50, unit):
    if unit == "millimolar":
        ic50 = math.pow(10, pic50)
    elif unit == "micromolar":
        ic50 = math.pow(10, 6-pic50)
    elif unit == "nanomolar":
        ic50 = math.pow(10, 9-pic50)
    return ic50

def ic50_to_pic50(ic50, unit):
    pic50 = None
    if unit == "millimolar":
        pic50 = -math.log10(ic50)
    elif unit == "micromolar":
        pic50 = 6 -math.log10(ic50 )
    elif unit == "nanomolar":
        pic50 = 9 - math.log10(ic50 )
    return pic50

def main():
    st.title("pIC50 to IC50 Converter")

    conversion_type = st.radio("Select conversion type:", ("pIC50 to IC50", "IC50 to pIC50"))

    default_unit = "nanometer"  # Default unit set to nanometer

    if conversion_type == "pIC50 to IC50":
        pic50 = st.number_input("Enter pIC50 value:", min_value=0.0)
        unit = st.selectbox("Select unit of concentration:", ("millimolar", "micromolar", "nanomolar"), index=2)
        if st.button("Convert"):
            ic50 = pic50_to_ic50(pic50, unit)
            st.success(f"IC50 value: {ic50:.6f} {unit}")

    elif conversion_type == "IC50 to pIC50":
        ic50 = st.number_input("Enter IC50 value:", min_value=0.0)
        unit = st.selectbox("Select unit of concentration:", ("millimolar", "micromolar", "nanomolar"), index=2)
        if st.button("Convert"):
            pic50 = ic50_to_pic50(ic50, unit)
            st.success(f"pIC50 value: {pic50:.2f}")

if __name__ == "__main__":
    main()
