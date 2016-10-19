function getRandomColour() {
    return colours[Math.floor(Math.random() * colours.length)];
}

var colours = ['Blue', 'Green', 'Purple', 'Red'];

$("#seqicon").attr('src', 'img/module_icons/Biojl_Seq_Icon_' + getRandomColour() + '.svg');
$("#alignicon").attr('src', 'img/module_icons/Biojl_Align_Icon_' + getRandomColour() + '.svg');
$("#intervalsicon").attr('src', 'img/module_icons/Biojl_Intervals_Icon_' + getRandomColour() + '.svg');
$("#structicon").attr('src', 'img/module_icons/Biojl_Struct_Icon_' + getRandomColour() + '.svg');
$("#varicon").attr('src', 'img/module_icons/Biojl_Var_Icon_' + getRandomColour() + '.svg');
$("#phyloicon").attr('src', 'img/module_icons/Biojl_Phylo_Icon_' + getRandomColour() + '.svg');
