#############################################################################
# Generated by PAGE version 4.7
# in conjunction with Tcl version 8.6
#    Oct 21, 2016 12:02:05 AM


set vTcl(actual_gui_bg) #d9d9d9
set vTcl(actual_gui_fg) #000000
set vTcl(actual_gui_menu_bg) #d9d9d9
set vTcl(actual_gui_menu_fg) #000000
set vTcl(complement_color) #d9d9d9
set vTcl(analog_color_p) #d9d9d9
set vTcl(analog_color_m) #d9d9d9
set vTcl(active_fg) #111111
#################################
#LIBRARY PROCEDURES
#


if {[info exists vTcl(sourcing)]} {

proc vTcl:project:info {} {
    set base .top36
    namespace eval ::widgets::$base {
        set dflt,origin 0
        set runvisible 1
    }
    namespace eval ::widgets_bindings {
        set tagslist _TopLevel
    }
    namespace eval ::vTcl::modules::main {
        set procs {
        }
        set compounds {
        }
        set projectType single
    }
}
}

#################################
# USER DEFINED PROCEDURES
#

#################################
# GENERATED GUI PROCEDURES
#

proc vTclWindow.top36 {base} {
    if {$base == ""} {
        set base .top36
    }
    if {[winfo exists $base]} {
        wm deiconify $base; return
    }
    set top $base
    ###################
    # CREATING WIDGETS
    ###################
    vTcl::widgets::core::toplevel::createCmd $top -class Toplevel \
        -background {#d9d9d9} -highlightbackground {#d9d9d9} \
        -highlightcolor black 
    wm focusmodel $top passive
    wm geometry $top 666x551+416+103
    update
    # set in toplevel.wgt.
    global vTcl
    set vTcl(save,dflt,origin) 0
    wm maxsize $top 3604 1185
    wm minsize $top 104 1
    wm overrideredirect $top 0
    wm resizable $top 1 1
    wm deiconify $top
    wm title $top "ZMCGUI"
    vTcl:DefineAlias "$top" "ZMCGUI" vTcl:Toplevel:WidgetProc "" 1
    button $top.but37 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -command OpenDataFile \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black -pady 0 \
        -text {input path to .mol2 file} 
    vTcl:DefineAlias "$top.but37" "Button1" vTcl:WidgetProc "ZMCGUI" 1
    text $top.tex39 \
        -background white -font TkTextFont -foreground black -height 36 \
        -highlightbackground {#d9d9d9} -highlightcolor black \
        -insertbackground black -selectbackground {#c4c4c4} \
        -selectforeground black -width 284 -wrap word 
    .top36.tex39 configure -font TkTextFont
    .top36.tex39 insert end text
    vTcl:DefineAlias "$top.tex39" "Text1" vTcl:WidgetProc "ZMCGUI" 1
    button $top.but40 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -command ImportModelforZMC \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black -pady 0 \
        -text {import model} 
    vTcl:DefineAlias "$top.but40" "Button2" vTcl:WidgetProc "ZMCGUI" 1
    label $top.lab41 \
        -activebackground {#f9f9f9} -activeforeground black \
        -background {#d9d9d9} -disabledforeground {#a3a3a3} \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -text {crystal size (a b c)} 
    vTcl:DefineAlias "$top.lab41" "Label1" vTcl:WidgetProc "ZMCGUI" 1
    entry $top.ent42 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black -textvariable 20 
    vTcl:DefineAlias "$top.ent42" "Entry1" vTcl:WidgetProc "ZMCGUI" 1
    entry $top.cpd44 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.cpd44" "Entry2" vTcl:WidgetProc "ZMCGUI" 1
    entry $top.cpd45 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.cpd45" "Entry3" vTcl:WidgetProc "ZMCGUI" 1
    label $top.lab48 \
        -activebackground {#f9f9f9} -activeforeground black \
        -background {#d9d9d9} -disabledforeground {#a3a3a3} \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -text {Unit cell parameters} 
    vTcl:DefineAlias "$top.lab48" "Label2" vTcl:WidgetProc "ZMCGUI" 1
    entry $top.ent49 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.ent49" "Entry4" vTcl:WidgetProc "ZMCGUI" 1
    entry $top.cpd50 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.cpd50" "Entry5" vTcl:WidgetProc "ZMCGUI" 1
    entry $top.cpd51 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.cpd51" "Entry6" vTcl:WidgetProc "ZMCGUI" 1
    entry $top.cpd52 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.cpd52" "Entry7" vTcl:WidgetProc "ZMCGUI" 1
    entry $top.cpd53 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.cpd53" "Entry8" vTcl:WidgetProc "ZMCGUI" 1
    entry $top.cpd54 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.cpd54" "Entry9" vTcl:WidgetProc "ZMCGUI" 1
    label $top.lab56 \
        -activebackground {#f9f9f9} -activeforeground black \
        -background {#d9d9d9} -disabledforeground {#a3a3a3} \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -text {loc zmat nmol} 
    vTcl:DefineAlias "$top.lab56" "Label3" vTcl:WidgetProc "ZMCGUI" 1
    entry $top.cpd58 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.cpd58" "Entry10" vTcl:WidgetProc "ZMCGUI" 1
    entry $top.cpd59 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.cpd59" "Entry11" vTcl:WidgetProc "ZMCGUI" 1
    entry $top.cpd60 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.cpd60" "Entry12" vTcl:WidgetProc "ZMCGUI" 1
    button $top.but61 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -command BuildRandomOcc \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black -pady 0 \
        -text {Build random occupancy file} 
    vTcl:DefineAlias "$top.but61" "Button3" vTcl:WidgetProc "ZMCGUI" 1
    button $top.but62 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -command RunZMATMaker \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black -pady 0 \
        -text {Run ZMAT Maker} 
    vTcl:DefineAlias "$top.but62" "Button4" vTcl:WidgetProc "ZMCGUI" 1
    button $top.but63 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -command ExitZMCGUI \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black -pady 0 \
        -text {Exit ZMCGUI} 
    vTcl:DefineAlias "$top.but63" "Button5" vTcl:WidgetProc "ZMCGUI" 1
    text $top.tex36 \
        -background white -font TkTextFont -foreground black -height 246 \
        -highlightbackground {#d9d9d9} -highlightcolor black \
        -insertbackground black -selectbackground {#c4c4c4} \
        -selectforeground black -width 304 -wrap word 
    .top36.tex36 configure -font TkTextFont
    .top36.tex36 insert end text
    entry $top.ent37 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable oxx.txt 
    vTcl:DefineAlias "$top.ent37" "Entry13" vTcl:WidgetProc "ZMCGUI" 1
    label $top.lab38 \
        -activebackground {#f9f9f9} -activeforeground black \
        -background {#d9d9d9} -disabledforeground {#a3a3a3} \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -text {occupancy filename } 
    vTcl:DefineAlias "$top.lab38" "Label4" vTcl:WidgetProc "ZMCGUI" 1
    button $top.but36 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -command GenQuickContacts \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black -pady 0 \
        -text {generate quick contacts} 
    vTcl:DefineAlias "$top.but36" "Button6" vTcl:WidgetProc "ZMCGUI" 1
    button $top.but38 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -command GenQuickZMCinput \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black -pady 0 \
        -text {generate quick ZMC input} 
    vTcl:DefineAlias "$top.but38" "Button7" vTcl:WidgetProc "ZMCGUI" 1
    button $top.but39 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -command ManageContacts \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black -pady 0 \
        -text {Manage Contacts} 
    vTcl:DefineAlias "$top.but39" "Button8" vTcl:WidgetProc "ZMCGUI" 1
    button $top.but41 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -command ZMCInputControl \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black -pady 0 \
        -text {ZMC input control} 
    vTcl:DefineAlias "$top.but41" "Button9" vTcl:WidgetProc "ZMCGUI" 1
    checkbutton $top.che42 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -disabledforeground {#a3a3a3} \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -justify left -text {output model for Diffuse} \
        -variable ZMC_option_Diffuse 
    vTcl:DefineAlias "$top.che42" "Checkbutton1" vTcl:WidgetProc "ZMCGUI" 1
    button $top.but44 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -command OpenDiffuseInputMenus \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black -pady 0 \
        -text {Open Diffuse input controls} 
    vTcl:DefineAlias "$top.but44" "Button10" vTcl:WidgetProc "ZMCGUI" 1
    button $top.but45 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -command RunZMC -disabledforeground {#a3a3a3} \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -pady 0 -text {Run ZMC} 
    vTcl:DefineAlias "$top.but45" "Button11" vTcl:WidgetProc "ZMCGUI" 1
    entry $top.ent36 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable Sites 
    vTcl:DefineAlias "$top.ent36" "Entry14" vTcl:WidgetProc "ZMCGUI" 1
    label $top.lab39 \
        -activebackground {#f9f9f9} -activeforeground black \
        -background {#d9d9d9} -disabledforeground {#a3a3a3} \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -text Sites 
    vTcl:DefineAlias "$top.lab39" "Label5" vTcl:WidgetProc "ZMCGUI" 1
    checkbutton $top.che36 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -disabledforeground {#a3a3a3} \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -justify left \
        -text {Auto relabel and add centriod} -variable DoRelabeling 
    vTcl:DefineAlias "$top.che36" "Checkbutton2" vTcl:WidgetProc "ZMCGUI" 1
    checkbutton $top.che37 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -disabledforeground {#a3a3a3} \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -justify left -text {Process without Centroid} \
        -variable NoCentroid 
    vTcl:DefineAlias "$top.che37" "Checkbutton3" vTcl:WidgetProc "ZMCGUI" 1
    button $top.but42 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -command ManageOccupancies \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black -pady 0 \
        -text {Manage Occupancies} 
    vTcl:DefineAlias "$top.but42" "Button12" vTcl:WidgetProc "ZMCGUI" 1
    label $top.lab36 \
        -activebackground {#f9f9f9} -activeforeground black \
        -background {#d9d9d9} -disabledforeground {#a3a3a3} \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -text {ZMCGUI v1.0 by Eric Chan (October 2016)} 
    vTcl:DefineAlias "$top.lab36" "Label6" vTcl:WidgetProc "ZMCGUI" 1
    button $top.but43 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -command WriteZMCRunScript \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black -pady 0 \
        -text {Write ZMC Run Script} 
    vTcl:DefineAlias "$top.but43" "Button13" vTcl:WidgetProc "ZMCGUI" 1
    checkbutton $top.che38 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -disabledforeground {#a3a3a3} \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -justify left -text {calc H} \
        -variable check_calcH 
    vTcl:DefineAlias "$top.che38" "Checkbutton4" vTcl:WidgetProc "ZMCGUI" 1
    checkbutton $top.che39 \
        -activebackground {#d9d9d9} -activeforeground {#000000} \
        -background {#d9d9d9} -disabledforeground {#a3a3a3} \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -justify left -text {calc C1} \
        -variable check_calcC1 
    vTcl:DefineAlias "$top.che39" "Checkbutton5" vTcl:WidgetProc "ZMCGUI" 1
    ###################
    # SETTING GEOMETRY
    ###################
    place $top.but37 \
        -in $top -x 30 -y 30 -width 141 -height 41 -anchor nw \
        -bordermode ignore 
    place $top.tex39 \
        -in $top -x 180 -y 30 -width 284 -height 36 -anchor nw \
        -bordermode ignore 
    place $top.but40 \
        -in $top -x 470 -y 30 -width 131 -height 31 -anchor nw \
        -bordermode ignore 
    place $top.lab41 \
        -in $top -x 60 -y 120 -width 101 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.ent42 \
        -in $top -x 180 -y 120 -width 44 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.cpd44 \
        -in $top -x 230 -y 120 -width 44 -height 19 -anchor nw \
        -bordermode inside 
    place $top.cpd45 \
        -in $top -x 280 -y 120 -width 44 -height 19 -anchor nw \
        -bordermode inside 
    place $top.lab48 \
        -in $top -x 50 -y 150 -anchor nw -bordermode ignore 
    place $top.ent49 \
        -in $top -x 170 -y 150 -width 74 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.cpd50 \
        -in $top -x 170 -y 180 -width 74 -height 19 -anchor nw \
        -bordermode inside 
    place $top.cpd51 \
        -in $top -x 170 -y 210 -width 74 -height 19 -anchor nw \
        -bordermode inside 
    place $top.cpd52 \
        -in $top -x 260 -y 150 -width 74 -height 19 -anchor nw \
        -bordermode inside 
    place $top.cpd53 \
        -in $top -x 260 -y 180 -width 74 -height 19 -anchor nw \
        -bordermode inside 
    place $top.cpd54 \
        -in $top -x 260 -y 210 -width 74 -height 19 -anchor nw \
        -bordermode inside 
    place $top.lab56 \
        -in $top -x 350 -y 120 -width 80 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.cpd58 \
        -in $top -x 430 -y 120 -width 44 -height 19 -anchor nw \
        -bordermode inside 
    place $top.cpd59 \
        -in $top -x 480 -y 120 -width 44 -height 19 -anchor nw \
        -bordermode inside 
    place $top.cpd60 \
        -in $top -x 530 -y 120 -width 44 -height 19 -anchor nw \
        -bordermode inside 
    place $top.but61 \
        -in $top -x 390 -y 230 -width 151 -height 31 -anchor nw \
        -bordermode ignore 
    place $top.but62 \
        -in $top -x 390 -y 190 -width 181 -height 31 -anchor nw \
        -bordermode ignore 
    place $top.but63 \
        -in $top -x 530 -y 480 -width 71 -height 51 -anchor nw \
        -bordermode ignore 
    place $top.tex36 \
        -in $top -x 30 -y 260 -width 304 -height 246 -anchor nw \
        -bordermode ignore 
    place $top.ent37 \
        -in $top -x 470 -y 150 -width 104 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.lab38 \
        -in $top -x 350 -y 150 -width 111 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.but36 \
        -in $top -x 390 -y 270 -width 131 -height 31 -anchor nw \
        -bordermode ignore 
    place $top.but38 \
        -in $top -x 390 -y 310 -width 131 -height 31 -anchor nw \
        -bordermode ignore 
    place $top.but39 \
        -in $top -x 530 -y 270 -width 101 -height 31 -anchor nw \
        -bordermode ignore 
    place $top.but41 \
        -in $top -x 530 -y 310 -width 101 -height 31 -anchor nw \
        -bordermode ignore 
    place $top.che42 \
        -in $top -x 380 -y 390 -width 167 -height 22 -anchor nw \
        -bordermode ignore 
    place $top.but44 \
        -in $top -x 390 -y 350 -width 161 -height 31 -anchor nw \
        -bordermode ignore 
    place $top.but45 \
        -in $top -x 390 -y 420 -width 91 -height 31 -anchor nw \
        -bordermode ignore 
    place $top.ent36 \
        -in $top -x 580 -y 120 -width 44 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.lab39 \
        -in $top -x 586 -y 96 -width 29 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.che36 \
        -in $top -x 280 -y 70 -width 177 -height 22 -anchor nw \
        -bordermode ignore 
    place $top.che37 \
        -in $top -x 274 -y 90 -width 167 -height 22 -anchor nw \
        -bordermode ignore 
    place $top.but42 \
        -in $top -x 550 -y 230 -width 111 -height 31 -anchor nw \
        -bordermode ignore 
    place $top.lab36 \
        -in $top -x 30 -y 7 -width 251 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.but43 \
        -in $top -x 490 -y 420 -width 121 -height 31 -anchor nw \
        -bordermode ignore 
    place $top.che38 \
        -in $top -x 573 -y 185 -width 57 -height 22 -anchor nw \
        -bordermode ignore 
    place $top.che39 \
        -in $top -x 575 -y 202 -width 57 -height 22 -anchor nw \
        -bordermode ignore 

    vTcl:FireEvent $base <<Ready>>
}

#############################################################################
## Binding tag:  _TopLevel

bind "_TopLevel" <<Create>> {
    if {![info exists _topcount]} {set _topcount 0}; incr _topcount
}
bind "_TopLevel" <<DeleteWindow>> {
    if {[set ::%W::_modal]} {
                vTcl:Toplevel:WidgetProc %W endmodal
            } else {
                destroy %W; if {$_topcount == 0} {exit}
            }
}
bind "_TopLevel" <Destroy> {
    if {[winfo toplevel %W] == "%W"} {incr _topcount -1}
}

Window show .
Window show .top36

