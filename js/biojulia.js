function getRandomIndex() {
    var randomIndex = Math.floor(Math.random() * bucket.length);
    return bucket.splice(randomIndex, 1)[0];
}

var bucket = [0, 1, 2, 3];
var ghIcons = ['Blue_GitHub.svg', 'Green_GitHub.svg', 'Purple_GitHub.svg', 'Red_GitHub.svg'];
var bkIcons = ['Blue_Book.svg', 'Green_Book.svg', 'Purple_Book.svg', 'Red_Book.svg'];
var gtIcons = ['Blue_Gitter.svg', 'Green_Gitter.svg', 'Purple_Gitter.svg', 'Red_Gitter.svg'];
var txtCols = ['#4266d5', '#3b972e', '#945bb0', '#c93d39'];

var ghidx = getRandomIndex();
$("#ghbutton").attr('src', 'img/github_icons/' + ghIcons[ghidx]);
$("#ghlink").css('color', txtCols[ghidx]);

var docidx = getRandomIndex();
$("#docbutton").attr('src', 'img/book_icons/' + bkIcons[docidx]);
$("#docslink").css('color', txtCols[docidx]);

var gttidx = getRandomIndex();
$("#gttbutton").attr('src', 'img/gitter_icons/' + gtIcons[gttidx]);
$("#gttlink").css('color', txtCols[gttidx]);
