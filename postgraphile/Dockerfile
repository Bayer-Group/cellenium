FROM node:alpine
LABEL description="Instant high-performance GraphQL API for your PostgreSQL database https://github.com/graphile/postgraphile"

# Install PostGraphile and PostGraphile connection filter plugin
RUN npm install -g postgraphile \
    npm install -g postgraphile-plugin-connection-filter \
    npm install -g @graphile-contrib/pg-simplify-inflector

EXPOSE 5000
ENTRYPOINT ["postgraphile", "-n", "0.0.0.0"]

