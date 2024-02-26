import streamlit as st
import math

def pic50_to_ic50(pic50):
    ic50 = 10 ** (-pic50)
    return ic50

def ic50_to_pic50(ic50):
    pic50 = -math.log10(ic50)
    return pic50

def main():
    st.title("🔢 PIC50 <-> IC50 Converter 🔢")
    st.write("This app converts between PIC50 and IC50 values.")

    conversion_mode = st.radio("Select conversion mode:", ("PIC50 to IC50", "IC50 to PIC50"))

    if conversion_mode == "PIC50 to IC50":
        pic50_input = st.number_input("Enter PIC50 value:", min_value=0.0, step=0.1, format="%.1f")
        if st.button("Convert"):
            ic50_result = pic50_to_ic50(pic50_input)
            st.write(f"The corresponding IC50 value is: {ic50_result:.5f}")
    elif conversion_mode == "IC50 to PIC50":
        ic50_input = st.number_input("Enter IC50 value:", min_value=0.0, step=0.0001, format="%.4f")
        if st.button("Convert"):
            pic50_result = ic50_to_pic50(ic50_input)
            st.write(f"The corresponding PIC50 value is: {pic50_result:.2f}")

if __name__ == "__main__":
    main()
