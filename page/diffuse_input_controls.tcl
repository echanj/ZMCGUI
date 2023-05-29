#############################################################################
# Generated by PAGE version 4.7
# in conjunction with Tcl version 8.6
#    Oct 22, 2016 12:15:58 PM


set vTcl(actual_gui_bg) wheat
set vTcl(actual_gui_fg) #000000
set vTcl(actual_gui_menu_bg) wheat
set vTcl(actual_gui_menu_fg) #000000
set vTcl(complement_color) #b2c9f4
set vTcl(analog_color_p) #eaf4b2
set vTcl(analog_color_m) #f4bcb2
set vTcl(active_fg) #111111
#################################
#LIBRARY PROCEDURES
#


if {[info exists vTcl(sourcing)]} {

proc vTcl:project:info {} {
    set base .top34
    namespace eval ::widgets::$base {
        set dflt,origin 0
        set runvisible 1
    }
    set site_3_0 $base.fra56
    set site_3_0 $base.fra62
    set site_3_0 $base.fra37
    set site_3_0 $base.fra39
    set site_3_0 $base.fra36
    set site_3_0 $base.fra38
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

proc vTclWindow.top34 {base} {
    if {$base == ""} {
        set base .top34
    }
    if {[winfo exists $base]} {
        wm deiconify $base; return
    }
    set top $base
    ###################
    # CREATING WIDGETS
    ###################
    vTcl::widgets::core::toplevel::createCmd $top -class Toplevel \
        -menu "$top.m40" -background {#8080c0} -highlightbackground wheat \
        -highlightcolor black 
    wm focusmodel $top passive
    wm geometry $top 600x445+502+145
    update
    # set in toplevel.wgt.
    global vTcl
    set vTcl(save,dflt,origin) 0
    wm maxsize $top 1905 1170
    wm minsize $top 104 1
    wm overrideredirect $top 0
    wm resizable $top 1 1
    wm deiconify $top
    wm title $top "Diffuse_Input_Controls"
    button $top.but34 \
        -activebackground {#f4bcb2} -activeforeground black \
        -background {#f5f5dedeb3b3} -command lambda:GenDiffuseInputFile() \
        -disabledforeground {#b8a786} -font font10 -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -pady 0 \
        -text {generate diffuse input} 
    vTcl:DefineAlias "$top.but34" "Button1" vTcl:WidgetProc "$top" 1
    button $top.but35 \
        -activebackground {#f4bcb2} -activeforeground black -background wheat \
        -command lambda:master.destroy() -disabledforeground {#b8a786} \
        -font font10 -foreground {#000000} -highlightbackground wheat \
        -highlightcolor black -pady 0 -text Close 
    vTcl:DefineAlias "$top.but35" "Button2" vTcl:WidgetProc "$top" 1
    label $top.lab36 \
        -activebackground {#ffffcd} -activeforeground black -anchor w \
        -background wheat -disabledforeground {#b8a786} -font font10 \
        -foreground {#000000} -highlightbackground wheat \
        -highlightcolor black -text Label -textvariable self.instance 
    vTcl:DefineAlias "$top.lab36" "Label1" vTcl:WidgetProc "$top" 1
    entry $top.ent36 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black -textvariable h1 
    vTcl:DefineAlias "$top.ent36" "Entry1" vTcl:WidgetProc "$top" 1
    entry $top.cpd37 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black -textvariable h3 
    vTcl:DefineAlias "$top.cpd37" "Entry2" vTcl:WidgetProc "$top" 1
    entry $top.cpd37.ent42 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black -textvariable h2 
    vTcl:DefineAlias "$top.cpd37.ent42" "Entry4" vTcl:WidgetProc "$top" 1
    entry $top.cpd38 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black -textvariable h3 
    vTcl:DefineAlias "$top.cpd38" "Entry3" vTcl:WidgetProc "$top" 1
    entry $top.cpd43 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black -textvariable k1 
    vTcl:DefineAlias "$top.cpd43" "Entry5" vTcl:WidgetProc "$top" 1
    entry $top.cpd44 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black -textvariable k2 
    vTcl:DefineAlias "$top.cpd44" "Entry6" vTcl:WidgetProc "$top" 1
    entry $top.cpd45 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black -textvariable k3 
    vTcl:DefineAlias "$top.cpd45" "Entry7" vTcl:WidgetProc "$top" 1
    entry $top.cpd46 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black -textvariable l1 
    vTcl:DefineAlias "$top.cpd46" "Entry8" vTcl:WidgetProc "$top" 1
    entry $top.cpd47 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black -textvariable l2 
    vTcl:DefineAlias "$top.cpd47" "Entry9" vTcl:WidgetProc "$top" 1
    entry $top.cpd48 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black -textvariable l3 
    vTcl:DefineAlias "$top.cpd48" "Entry10" vTcl:WidgetProc "$top" 1
    label $top.lab37 \
        -activebackground {#f9f9f9} -activeforeground black -background wheat \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black -text Origin 
    vTcl:DefineAlias "$top.lab37" "Label2" vTcl:WidgetProc "$top" 1
    label $top.cpd39 \
        -activebackground {#f9f9f9} -activeforeground black -background wheat \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black \
        -text {Vert. axis} 
    vTcl:DefineAlias "$top.cpd39" "Label3" vTcl:WidgetProc "$top" 1
    menu $top.m40 \
        -activebackground {#f4bcb2} -activeforeground black -background wheat \
        -font TkMenuFont -foreground {#000000} -tearoff 0 
    label $top.cpd41 \
        -activebackground {#f9f9f9} -activeforeground black -background wheat \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black \
        -text {Horz. axis} 
    vTcl:DefineAlias "$top.cpd41" "Label4" vTcl:WidgetProc "$top" 1
    button $top.but42 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -command run_diffuse_job \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -pady 0 \
        -text {Run diffuse Job} 
    vTcl:DefineAlias "$top.but42" "Button3" vTcl:WidgetProc "$top" 1
    entry $top.ent43 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable dsname 
    vTcl:DefineAlias "$top.ent43" "Entry11" vTcl:WidgetProc "$top" 1
    label $top.lab44 \
        -activebackground {#f9f9f9} -activeforeground black -background wheat \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black \
        -text {name (e.g. h0l)} 
    vTcl:DefineAlias "$top.lab44" "Label5" vTcl:WidgetProc "$top" 1
    entry $top.ent45 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black -textvariable Lam 
    vTcl:DefineAlias "$top.ent45" "Entry12" vTcl:WidgetProc "$top" 1
    entry $top.ent46 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.ent46" "Entry13" vTcl:WidgetProc "$top" 1
    label $top.lab47 \
        -activebackground {#f9f9f9} -activeforeground black -background wheat \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black -text Lam,Thmax 
    vTcl:DefineAlias "$top.lab47" "Label6" vTcl:WidgetProc "$top" 1
    text $top.tex48 \
        -background white -font TkTextFont -foreground black -height 66 \
        -highlightbackground wheat -highlightcolor black \
        -insertbackground black -selectbackground {#c4c4c4} \
        -selectforeground black -width 344 -wrap word 
    .top34.tex48 configure -font TkTextFont
    .top34.tex48 insert end text
    vTcl:DefineAlias "$top.tex48" "Text1" vTcl:WidgetProc "$top" 1
    entry $top.ent37 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable Lot_a 
    vTcl:DefineAlias "$top.ent37" "Entry14" vTcl:WidgetProc "$top" 1
    label $top.lab38 \
        -activebackground {#f9f9f9} -activeforeground black -background wheat \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black \
        -text {Lotsizes a b c} 
    vTcl:DefineAlias "$top.lab38" "Label7" vTcl:WidgetProc "$top" 1
    entry $top.cpd40 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.cpd40" "Entry15" vTcl:WidgetProc "$top" 1
    entry $top.cpd42 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable lotc 
    vTcl:DefineAlias "$top.cpd42" "Entry16" vTcl:WidgetProc "$top" 1
    entry $top.ent44 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.ent44" "Entry17" vTcl:WidgetProc "$top" 1
    label $top.lab45 \
        -activebackground {#f9f9f9} -activeforeground black -background wheat \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black \
        -text {num. lots } 
    vTcl:DefineAlias "$top.lab45" "Label8" vTcl:WidgetProc "$top" 1
    entry $top.ent47 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.ent47" "Entry18" vTcl:WidgetProc "$top" 1
    label $top.lab48 \
        -activebackground {#f9f9f9} -activeforeground black -background wheat \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black \
        -text {atom sites} 
    vTcl:DefineAlias "$top.lab48" "Label9" vTcl:WidgetProc "$top" 1
    entry $top.ent50 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black 
    vTcl:DefineAlias "$top.ent50" "Entry19" vTcl:WidgetProc "$top" 1
    label $top.lab51 \
        -activebackground {#f9f9f9} -activeforeground black -background wheat \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black \
        -text {atom types} 
    vTcl:DefineAlias "$top.lab51" "Label10" vTcl:WidgetProc "$top" 1
    checkbutton $top.che55 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text {Do PBC ?} -variable check_doPBC 
    vTcl:DefineAlias "$top.che55" "Checkbutton1" vTcl:WidgetProc "$top" 1
    frame $top.fra56 \
        -borderwidth 2 -relief groove -background wheat -height 85 \
        -highlightbackground {#d9d9d9} -highlightcolor black -width 145 
    vTcl:DefineAlias "$top.fra56" "Frame1" vTcl:WidgetProc "$top" 1
    set site_3_0 $top.fra56
    radiobutton $site_3_0.rad57 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text None -value 1 -variable {} 
    vTcl:DefineAlias "$site_3_0.rad57" "Radiobutton1" vTcl:WidgetProc "$top" 1
    radiobutton $site_3_0.rad58 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text aver -variable Bragg_treatment 
    vTcl:DefineAlias "$site_3_0.rad58" "Radiobutton2" vTcl:WidgetProc "$top" 1
    radiobutton $site_3_0.rad59 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text 5% -value 0 -variable Bragg_treatment 
    vTcl:DefineAlias "$site_3_0.rad59" "Radiobutton3" vTcl:WidgetProc "$top" 1
    label $site_3_0.lab60 \
        -activebackground {#f9f9f9} -activeforeground black -background wheat \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black \
        -text {Bragg removal
options} 
    vTcl:DefineAlias "$site_3_0.lab60" "Label11" vTcl:WidgetProc "$top" 1
    radiobutton $site_3_0.rad61 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text All -value 0 -variable Bragg_treatment 
    vTcl:DefineAlias "$site_3_0.rad61" "Radiobutton4" vTcl:WidgetProc "$top" 1
    place $site_3_0.rad57 \
        -in $site_3_0 -x 10 -y 30 -width 55 -height 32 -anchor nw \
        -bordermode ignore 
    place $site_3_0.rad58 \
        -in $site_3_0 -x 69 -y 36 -width 55 -height 22 -anchor nw \
        -bordermode ignore 
    place $site_3_0.rad59 \
        -in $site_3_0 -x 2 -y 51 -width 65 -height 32 -anchor nw \
        -bordermode ignore 
    place $site_3_0.lab60 \
        -in $site_3_0 -x 7 -y 4 -width 121 -height 29 -anchor nw \
        -bordermode ignore 
    place $site_3_0.rad61 \
        -in $site_3_0 -x 64 -y 56 -width 55 -height 22 -anchor nw \
        -bordermode ignore 
    frame $top.fra62 \
        -borderwidth 2 -relief groove -background wheat -height 125 \
        -highlightbackground {#d9d9d9} -highlightcolor black -width 125 
    vTcl:DefineAlias "$top.fra62" "Frame2" vTcl:WidgetProc "$top" 1
    set site_3_0 $top.fra62
    label $site_3_0.lab63 \
        -activebackground {#f9f9f9} -activeforeground black -background wheat \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black \
        -text {Which Atoms ?} 
    vTcl:DefineAlias "$site_3_0.lab63" "Label12" vTcl:WidgetProc "$top" 1
    checkbutton $site_3_0.che64 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text Carbon -variable checkCarbon 
    vTcl:DefineAlias "$site_3_0.che64" "Checkbutton2" vTcl:WidgetProc "$top" 1
    checkbutton $site_3_0.che65 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text Oxygen -variable checkOxygen 
    vTcl:DefineAlias "$site_3_0.che65" "Checkbutton3" vTcl:WidgetProc "$top" 1
    checkbutton $site_3_0.che66 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text Nitrogen -variable checkNitrogen 
    vTcl:DefineAlias "$site_3_0.che66" "Checkbutton4" vTcl:WidgetProc "$top" 1
    checkbutton $site_3_0.che67 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text Sulfur -variable checkSulfur 
    vTcl:DefineAlias "$site_3_0.che67" "Checkbutton5" vTcl:WidgetProc "$top" 1
    checkbutton $site_3_0.che68 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text Chlorine -variable checkChlorine 
    vTcl:DefineAlias "$site_3_0.che68" "Checkbutton6" vTcl:WidgetProc "$top" 1
    checkbutton $site_3_0.che36 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text Fluorine -variable check_Fluorine 
    vTcl:DefineAlias "$site_3_0.che36" "Checkbutton11" vTcl:WidgetProc "$top" 1
    place $site_3_0.lab63 \
        -in $site_3_0 -x 20 -y 5 -width 81 -height 19 -anchor nw \
        -bordermode ignore 
    place $site_3_0.che64 \
        -in $site_3_0 -x 12 -y 22 -width 57 -height 22 -anchor nw \
        -bordermode ignore 
    place $site_3_0.che65 \
        -in $site_3_0 -x 13 -y 38 -width 57 -height 22 -anchor nw \
        -bordermode ignore 
    place $site_3_0.che66 \
        -in $site_3_0 -x 15 -y 56 -width 57 -height 22 -anchor nw \
        -bordermode ignore 
    place $site_3_0.che67 \
        -in $site_3_0 -x 9 -y 74 -width 57 -height 22 -anchor nw \
        -bordermode ignore 
    place $site_3_0.che68 \
        -in $site_3_0 -x 9 -y 91 -width 67 -height 22 -anchor nw \
        -bordermode ignore 
    place $site_3_0.che36 \
        -in $site_3_0 -x 13 -y 107 -width 57 -height 22 -anchor nw \
        -bordermode ignore 
    frame $top.fra37 \
        -borderwidth 2 -relief groove -background wheat -height 55 \
        -highlightbackground {#d9d9d9} -highlightcolor black -width 285 
    vTcl:DefineAlias "$top.fra37" "Frame3" vTcl:WidgetProc "$top" 1
    set site_3_0 $top.fra37
    checkbutton $site_3_0.che38 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text Twofold -variable check_twofold 
    vTcl:DefineAlias "$site_3_0.che38" "Checkbutton7" vTcl:WidgetProc "$top" 1
    checkbutton $site_3_0.che39 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text {horizontal
mirrior} -variable check_horiz_mirror 
    vTcl:DefineAlias "$site_3_0.che39" "Checkbutton8" vTcl:WidgetProc "$top" 1
    checkbutton $site_3_0.che40 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text {vertical
mirror} -variable check_vert_mirror 
    vTcl:DefineAlias "$site_3_0.che40" "Checkbutton9" vTcl:WidgetProc "$top" 1
    checkbutton $site_3_0.che41 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text Normalize -variable check_normalize 
    vTcl:DefineAlias "$site_3_0.che41" "Checkbutton10" vTcl:WidgetProc "$top" 1
    entry $site_3_0.ent43 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable norm_val 
    vTcl:DefineAlias "$site_3_0.ent43" "Entry20" vTcl:WidgetProc "$top" 1
    label $site_3_0.lab44 \
        -activebackground {#f9f9f9} -activeforeground black -background wheat \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black \
        -text {bin2gray options } 
    vTcl:DefineAlias "$site_3_0.lab44" "Label13" vTcl:WidgetProc "$top" 1
    place $site_3_0.che38 \
        -in $site_3_0 -x 6 -y 23 -width 66 -height 22 -anchor nw \
        -bordermode ignore 
    place $site_3_0.che39 \
        -in $site_3_0 -x 71 -y 16 -width 75 -height 35 -anchor nw \
        -bordermode ignore 
    place $site_3_0.che40 \
        -in $site_3_0 -x 142 -y 15 -width 63 -height 35 -anchor nw \
        -bordermode ignore 
    place $site_3_0.che41 \
        -in $site_3_0 -x 199 -y 5 -width 74 -height 22 -anchor nw \
        -bordermode ignore 
    place $site_3_0.ent43 \
        -in $site_3_0 -x 207 -y 28 -width 64 -height 19 -anchor nw \
        -bordermode ignore 
    place $site_3_0.lab44 \
        -in $site_3_0 -x 7 -y 4 -width 101 -height 19 -anchor nw \
        -bordermode ignore 
    button $top.but36 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -command GetLotSizes -disabledforeground {#a3a3a3} \
        -foreground {#000000} -highlightbackground wheat \
        -highlightcolor black -pady 0 -text {Get Lot sizes} 
    vTcl:DefineAlias "$top.but36" "Button4" vTcl:WidgetProc "$top" 1
    entry $top.ent38 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable GetLotSize_iter 
    vTcl:DefineAlias "$top.ent38" "Entry21" vTcl:WidgetProc "$top" 1
    frame $top.fra39 \
        -borderwidth 2 -relief groove -background wheat -height 35 \
        -highlightbackground {#d9d9d9} -highlightcolor black -width 195 
    vTcl:DefineAlias "$top.fra39" "Frame4" vTcl:WidgetProc "$top" 1
    set site_3_0 $top.fra39
    radiobutton $site_3_0.rad40 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text h0l -value 0 -variable radio_hkl 
    vTcl:DefineAlias "$site_3_0.rad40" "Radiobutton5" vTcl:WidgetProc "$top" 1
    radiobutton $site_3_0.rad41 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text hk0 -value 1 -variable radio_hkl 
    vTcl:DefineAlias "$site_3_0.rad41" "Radiobutton6" vTcl:WidgetProc "$top" 1
    radiobutton $site_3_0.rad42 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text 0kl -value 2 -variable radio_hkl 
    vTcl:DefineAlias "$site_3_0.rad42" "Radiobutton7" vTcl:WidgetProc "$top" 1
    place $site_3_0.rad40 \
        -in $site_3_0 -x 5 -y 6 -width 55 -height 22 -anchor nw \
        -bordermode ignore 
    place $site_3_0.rad41 \
        -in $site_3_0 -x 64 -y 7 -width 55 -height 22 -anchor nw \
        -bordermode ignore 
    place $site_3_0.rad42 \
        -in $site_3_0 -x 128 -y 6 -width 55 -height 22 -anchor nw \
        -bordermode ignore 
    button $top.but44 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -command makeDiffuseRunScript \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -pady 0 \
        -text {make diffuse run script} 
    vTcl:DefineAlias "$top.but44" "Button5" vTcl:WidgetProc "$top" 1
    frame $top.fra36 \
        -borderwidth 2 -relief groove -background wheat -height 35 \
        -highlightbackground {#d9d9d9} -highlightcolor black -width 195 
    vTcl:DefineAlias "$top.fra36" "Frame5" vTcl:WidgetProc "$top" 1
    set site_3_0 $top.fra36
    entry $site_3_0.ent37 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable PixelsX 
    vTcl:DefineAlias "$site_3_0.ent37" "Entry22" vTcl:WidgetProc "$top" 1
    entry $site_3_0.ent38 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable PixelsY 
    vTcl:DefineAlias "$site_3_0.ent38" "Entry23" vTcl:WidgetProc "$top" 1
    entry $site_3_0.ent39 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -highlightbackground {#d9d9d9} \
        -highlightcolor black -insertbackground black \
        -selectbackground {#c4c4c4} -selectforeground black \
        -textvariable PixelsZ 
    vTcl:DefineAlias "$site_3_0.ent39" "Entry24" vTcl:WidgetProc "$top" 1
    label $site_3_0.lab40 \
        -activebackground {#f9f9f9} -activeforeground black -background wheat \
        -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground {#d9d9d9} -highlightcolor black \
        -text {Pixels (XYZ)} 
    vTcl:DefineAlias "$site_3_0.lab40" "Label14" vTcl:WidgetProc "$top" 1
    place $site_3_0.ent37 \
        -in $site_3_0 -x 7 -y 8 -width 34 -height 19 -anchor nw \
        -bordermode ignore 
    place $site_3_0.ent38 \
        -in $site_3_0 -x 47 -y 8 -width 34 -height 19 -anchor nw \
        -bordermode ignore 
    place $site_3_0.ent39 \
        -in $site_3_0 -x 85 -y 8 -width 34 -height 19 -anchor nw \
        -bordermode ignore 
    place $site_3_0.lab40 \
        -in $site_3_0 -x 124 -y 8 -width 62 -height 19 -anchor nw \
        -bordermode ignore 
    frame $top.fra38 \
        -borderwidth 2 -relief groove -background wheat -height 55 -width 195 
    vTcl:DefineAlias "$top.fra38" "Frame6" vTcl:WidgetProc "$top" 1
    set site_3_0 $top.fra38
    checkbutton $site_3_0.che39 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -highlightbackground wheat -highlightcolor black -justify left \
        -text {Apply cutoffs
for viewing} -variable check_vMinMax 
    vTcl:DefineAlias "$site_3_0.che39" "Checkbutton12" vTcl:WidgetProc "$top" 1
    entry $site_3_0.ent40 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -insertbackground black -textvariable pgm_vMin 
    vTcl:DefineAlias "$site_3_0.ent40" "Entry25" vTcl:WidgetProc "$top" 1
    entry $site_3_0.ent41 \
        -background white -disabledforeground {#a3a3a3} -font TkFixedFont \
        -foreground {#000000} -insertbackground black -textvariable pgm_vMax 
    vTcl:DefineAlias "$site_3_0.ent41" "Entry26" vTcl:WidgetProc "$top" 1
    label $site_3_0.lab42 \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -text Min 
    vTcl:DefineAlias "$site_3_0.lab42" "Label15" vTcl:WidgetProc "$top" 1
    label $site_3_0.lab43 \
        -background wheat -disabledforeground {#a3a3a3} -foreground {#000000} \
        -text Max 
    vTcl:DefineAlias "$site_3_0.lab43" "Label16" vTcl:WidgetProc "$top" 1
    place $site_3_0.che39 \
        -in $site_3_0 -x 4 -y 6 -width 87 -height 42 -anchor nw \
        -bordermode ignore 
    place $site_3_0.ent40 \
        -in $site_3_0 -x 97 -y 5 -width 64 -height 19 -anchor nw \
        -bordermode ignore 
    place $site_3_0.ent41 \
        -in $site_3_0 -x 97 -y 30 -width 64 -height 19 -anchor nw \
        -bordermode ignore 
    place $site_3_0.lab42 \
        -in $site_3_0 -x 165 -y 5 -width 22 -height 19 -anchor nw \
        -bordermode ignore 
    place $site_3_0.lab43 \
        -in $site_3_0 -x 164 -y 29 -width 26 -height 19 -anchor nw \
        -bordermode ignore 
    button $top.but45 \
        -activebackground {#f4bcb2} -activeforeground {#000000} \
        -background wheat -command ViewPGMImage -disabledforeground {#a3a3a3} \
        -foreground {#000000} -highlightbackground wheat \
        -highlightcolor black -pady 0 -text {View Image} 
    vTcl:DefineAlias "$top.but45" "Button6" vTcl:WidgetProc "$top" 1
    ###################
    # SETTING GEOMETRY
    ###################
    place $top.but34 \
        -in $top -x 380 -y 350 -width 168 -height 30 -anchor nw \
        -bordermode ignore 
    place $top.but35 \
        -in $top -x 540 -y 410 -width 55 -height 30 -anchor nw \
        -bordermode ignore 
    place $top.lab36 \
        -in $top -x 160 -y 20 -width 283 -height 27 -anchor nw \
        -bordermode ignore 
    place $top.ent36 \
        -in $top -x 14 -y 101 -width 34 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.cpd37 \
        -in $top -x 56 -y 101 -width 34 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.cpd37.ent42 \
        -in $top.cpd37 -x 0 -y 0 -anchor nw -bordermode ignore 
    place $top.cpd38 \
        -in $top -x 98 -y 101 -width 34 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.cpd43 \
        -in $top -x 14 -y 124 -width 34 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.cpd44 \
        -in $top -x 56 -y 124 -width 34 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.cpd45 \
        -in $top -x 98 -y 124 -width 34 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.cpd46 \
        -in $top -x 14 -y 146 -width 34 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.cpd47 \
        -in $top -x 56 -y 145 -width 34 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.cpd48 \
        -in $top -x 97 -y 146 -width 34 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.lab37 \
        -in $top -x 140 -y 100 -width 61 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.cpd39 \
        -in $top -x 139 -y 124 -width 61 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.cpd41 \
        -in $top -x 139 -y 147 -width 61 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.but42 \
        -in $top -x 380 -y 384 -width 91 -height 31 -anchor nw \
        -bordermode ignore 
    place $top.ent43 \
        -in $top -x 18 -y 70 -width 84 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.lab44 \
        -in $top -x 107 -y 70 -width 91 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.ent45 \
        -in $top -x 14 -y 171 -width 54 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.ent46 \
        -in $top -x 76 -y 171 -width 54 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.lab47 \
        -in $top -x 139 -y 170 -width 61 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.tex48 \
        -in $top -x 24 -y 351 -width 344 -height 66 -anchor nw \
        -bordermode ignore 
    place $top.ent37 \
        -in $top -x 225 -y 71 -width 54 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.lab38 \
        -in $top -x 413 -y 71 -width 81 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.cpd40 \
        -in $top -x 286 -y 71 -width 54 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.cpd42 \
        -in $top -x 345 -y 71 -width 54 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.ent44 \
        -in $top -x 225 -y 97 -width 54 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.lab45 \
        -in $top -x 287 -y 97 -width 71 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.ent47 \
        -in $top -x 226 -y 121 -width 54 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.lab48 \
        -in $top -x 286 -y 122 -width 71 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.ent50 \
        -in $top -x 226 -y 146 -width 54 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.lab51 \
        -in $top -x 286 -y 146 -width 71 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.che55 \
        -in $top -x 227 -y 171 -width 77 -height 22 -anchor nw \
        -bordermode ignore 
    place $top.fra56 \
        -in $top -x 225 -y 196 -width 145 -height 85 -anchor nw \
        -bordermode ignore 
    place $top.fra62 \
        -in $top -x 376 -y 156 -width 125 -height 125 -anchor nw \
        -bordermode ignore 
    place $top.fra37 \
        -in $top -x 226 -y 286 -width 285 -height 55 -anchor nw \
        -bordermode ignore 
    place $top.but36 \
        -in $top -x 414 -y 95 -width 71 -height 21 -anchor nw \
        -bordermode ignore 
    place $top.ent38 \
        -in $top -x 489 -y 95 -width 34 -height 19 -anchor nw \
        -bordermode ignore 
    place $top.fra39 \
        -in $top -x 10 -y 200 -width 195 -height 35 -anchor nw \
        -bordermode ignore 
    place $top.but44 \
        -in $top -x 380 -y 420 -width 151 -height 21 -anchor nw \
        -bordermode ignore 
    place $top.fra36 \
        -in $top -x 10 -y 240 -width 195 -height 35 -anchor nw \
        -bordermode ignore 
    place $top.fra38 \
        -in $top -x 10 -y 280 -width 195 -height 55 -anchor nw \
        -bordermode ignore 
    place $top.but45 \
        -in $top -x 474 -y 384 -width 61 -height 31 -anchor nw \
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
Window show .top34

