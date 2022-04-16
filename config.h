#define PLUG_NAME "Tablitsa2"
#define PLUG_MFR "RadioFarmer"
#define PLUG_VERSION_HEX 0x02000000
#define PLUG_VERSION_STR "2.0.0.0"
#define PLUG_UNIQUE_ID 'U8jr'
#define PLUG_MFR_ID 'RdFm'
#define PLUG_URL_STR "https://github.com/radiofarmer/Tablitsa2-Synthesizer"
#define PLUG_EMAIL_STR ""
#define PLUG_COPYRIGHT_STR "Copyright Bram Osterhout 2021, MIT License"
#define PLUG_CLASS_NAME Tablitsa2

#define BUNDLE_NAME "Tablitsa2 (Beta)"
#define BUNDLE_MFR "RadioFarmer"
#define BUNDLE_DOMAIN "com"

#define PLUG_CHANNEL_IO "0-2"
#define SHARED_RESOURCES_SUBPATH "Tablitsa2"

#define PLUG_LATENCY 0
#define PLUG_TYPE 1
#define PLUG_DOES_MIDI_IN 1
#define PLUG_DOES_MIDI_OUT 1
#define PLUG_DOES_MPE 1
#define PLUG_DOES_STATE_CHUNKS 1
#define PLUG_HAS_UI 1
#define PLUG_WIDTH 1400
#define PLUG_HEIGHT 770
#define PLUG_FPS 60
#define PLUG_SHARED_RESOURCES 0
#define PLUG_HOST_RESIZE 1

#define AUV2_ENTRY Tablitsa2_Entry
#define AUV2_ENTRY_STR "Tablitsa2_Entry"
#define AUV2_FACTORY Tablitsa2_Factory
#define AUV2_VIEW_CLASS Tablitsa2_View
#define AUV2_VIEW_CLASS_STR "Tablitsa2_View"

#define AAX_TYPE_IDS 'IPI1', 'IPI2'
#define AAX_PLUG_MFR_STR "Acme"
#define AAX_PLUG_NAME_STR "Tablitsa2\nIPIS"
#define AAX_DOES_AUDIOSUITE 0
#define AAX_PLUG_CATEGORY_STR "Synth"

#define VST3_SUBCATEGORY "Instrument|Synth"

#define APP_NUM_CHANNELS 2
#define APP_N_VECTOR_WAIT 0
#define APP_MULT 1
#define APP_COPY_AUV3 0
#define APP_SIGNAL_VECTOR_SIZE 64

// NB: Remember to edit main.rc compile-time directives
#define ROBOTO_FN "Roboto-Regular.ttf"
#define BG_FN "tablitsa_background.png"
#define PERIODIC_TABLE_FN "periodic_table.svg"
#define LOGO_FN "logo.png"
#define ROUTING_FN "routing_diagram.svg"
