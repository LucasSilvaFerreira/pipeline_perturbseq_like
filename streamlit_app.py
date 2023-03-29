from st_on_hover_tabs import on_hover_tabs
import streamlit as st
from streamlit_extras.metric_cards import style_metric_cards



style_metric_cards()

st.header("Single cell Perturb-seq like REPORT")
st.markdown('<style>' + open('./style.css').read() + '</style>', unsafe_allow_html=True)

col1, col2 = st.columns(2)

with st.sidebar:
    tabs = on_hover_tabs(tabName=['page1', 'page2', 'page3'], 
                         iconName=['dashboard', 'money', 'economy'], default_choice=0)

if tabs =='page1':
    st.title("Sample name")
    st.write('Name of option is {}'.format(tabs))

    with col1:

      st.header(":blue[Guide Mapping]")
      st.metric(label="Total number of cells", value=1000, delta="number of barcodes:100")


      #st.subheader(':blue[1000]')
      v1a, v1b = st.columns(2)
      with v1a:
        st.write('Total number of reads')
        st.write('Mapped Reads') 
        st.write('% of mapped reads')

      with v1b:
        st.write('1000')
        st.write('1000')
        st.write('%100')

      st.success('This is a success message!', icon="✅")

      st.subheader('Total number of reads')
      st.image("https://www150.statcan.gc.ca/edu/power-pouvoir/c-g/c-g05-2-1-eng.png")
      st.subheader('Total number of reads')
      st.image("https://miro.medium.com/v2/resize:fit:841/1*-MSLBij5vknV5Tm-G0Rxfw.png")
      st.subheader('Total number of reads')
      st.image("https://miro.medium.com/v2/resize:fit:841/1*-MSLBij5vknV5Tm-G0Rxfw.png")


    with col2:
      st.header(":red[scRNAseq Mapping]")
      st.metric(label="Total number of cells", value=1000, delta="number of barcodes:100")

      st.write('Total number of reads: 1000')
      st.write('Mapped Reads: 1000')
      st.write('% of mapped reads: 89%')
      st.success('This is a success message!', icon="✅")


      #st.image("https://static.streamlit.io/examples/dog.jpg")


    

elif tabs == 'page2':
    st.title("page2")
    st.write('Name of option is {}'.format(tabs))

elif tabs == 'page3':
    st.title("page3")
    st.write('Name of option is {}'.format(tabs))
    
    
