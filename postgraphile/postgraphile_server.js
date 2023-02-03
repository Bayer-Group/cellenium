// see https://www.graphile.org/postgraphile/usage-library/
const express = require("express");
const {postgraphile} = require("postgraphile");
const app = express();
var set = require('lodash.set');

// see express-jwt implementation example
var inspectJwtMiddleware = function (req, res, next) {
    // console.log('all headers', req.headers);
    if (process.env.AUTH_HEADER_NAME) {
        header_value = req.headers[process.env.AUTH_HEADER_NAME];
        if (!header_value && process.env.AUTH_HEADER_MOCK_VALUE) {
            header_value = process.env.AUTH_HEADER_MOCK_VALUE;
        }
        if (header_value) {
            set(req, 'auth_header_value', header_value);
        }
    }
    next();
}

// Apply checkJwt to our graphql endpoint
app.use("/postgraphile/", inspectJwtMiddleware);

app.use(
    postgraphile(
        process.env.DATABASE_URL,
        "public",
        {
            pgSettings: req => {
                const settings = {};
                if (req.auth_header_value) {
                    settings["postgraphile.auth_header_value"] = req.auth_header_value;
                }
                // console.log('settings', settings);
                return settings;
            },
            watchPg: !!process.env.OWNER_DATABASE_URL,
            ownerConnectionString: process.env.OWNER_DATABASE_URL,
            retryOnInitFail: true,
            graphqlRoute: "/postgraphile/",
            graphiql: true,
            graphiqlRoute: "/postgraphile/graphiql",
            enhanceGraphiql: true,
            allowExplain: true,
            dynamicJson: true,
            showErrorStack: "json",
            extendedErrors: ["hint", "detail", "errcode"],
            simpleCollections: "only",
            enableQueryBatching: true,
            appendPlugins: [require("@graphile-contrib/pg-simplify-inflector"),
                require("postgraphile-plugin-connection-filter")]
        }
    )
);

app.listen(process.env.PORT || 5000);